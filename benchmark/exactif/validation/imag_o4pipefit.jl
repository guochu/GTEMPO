# =============================================================================
# imag_o4pipefit.jl —— 表示精确的 I 链求解：真实管线循环内牛顿
#
# 原理（表示依赖发现后）：复合量依赖 MPS 因子化，多项式级拟合原理性失效。
# 本脚本把 I 链 = ∏_a exp(GTerm(row_a))（7×7 对族、行序）的系数 Δ̃ 直接在
# **真实管线**里牛顿求解：残差 = Z/gf 复合 vs Fock 锚。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4pipefit.jl
# =============================================================================
using Printf
using LinearAlgebra
using Random
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const NL = 2
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

annih_op_f(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end
function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op_f(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) || (Hq .+= h[i, j] .* (ops[i]' * ops[j]))
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Z = Matrix{ComplexF64}(I, dim, dim); Zn = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z; Zn = W * Zn
    end
    Zn = ops[1]' * ops[1] * Zn
    return tr(Z), tr(Zn) / tr(Z)
end

const ORDER = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
const SCHED = NativeExact.native_schedule(ORDER, NL, BETA / NL)
const RLINKS = vcat(SCHED.links[1] + SCHED.links[end], SCHED.links[2:end-1])
const NL_INT = length(RLINKS)
const CELLS = SCHED.cells
const LAT = GTEMPO.GrassmannLattice(N=NL_INT, δτ=BETA / NL_INT, contour=:imag, bands=1)
const K = LAT.k
const SL = collect(1:K)
const BARP = [LAT[GrassmannContourIndex(s, :τ, 1, true)] for s in SL]
const ETAP = [LAT[GrassmannContourIndex(s, :τ, 1, false)] for s in SL]

"""K 链（交错：门@切片 g+1，分辨率 −1，BC）。"""
function build_kchain(eps_d)
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for g in 1:NL_INT
        iszero(eps_d) || begin
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(BARP[g+1 <= K ? g+1 : g], ETAP[g+1 <= K ? g+1 : g];
                                 coeff=ComplexF64(exp(-RLINKS[g] * eps_d) - 1.0))), mps)
        end
        s1 = g + 1
        s2 = g == NL_INT ? 2 : g + 2
        rb = g == NL_INT ? BARP[2] : BARP[s1]
        re = g == NL_INT ? ETAP[2] : ETAP[s2]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(rb, re; coeff=ComplexF64(-1.0))), mps)
    end
    mps = GTEMPO.apply!(exp(GTEMPO.GTerm(BARP[1], ETAP[K]; coeff=ComplexF64(-1.0))), mps)
    return mps
end

"""I-MPS：7×7 Δ̃ 行序注入。"""
function build_imps(x)
    m = GTEMPO.vacuumstate(ComplexF64, LAT)
    for a in 1:K
        for b in 1:K
            re = x[2*(K*(a-1)+b)-1]; im = x[2*(K*(a-1)+b)]
            abs(re) < 1e-14 && abs(im) < 1e-14 && continue
            m = GTEMPO.apply!(exp(GTEMPO.GTerm(BARP[a], ETAP[b];
                                 coeff=ComplexF64(re, im))), m)
        end
    end
    return m
end

"""复合：Z 与 1−g_end（真实 integrate/gf）。"""
function composite(x, eps_d; mpsK_cache=nothing)
    mpsK = mpsK_cache === nothing ? build_kchain(eps_d) : mpsK_cache
    mpsI = build_imps(x)
    z = GTEMPO.integrate(LAT, mpsK, mpsI) * x[end]
    b1 = ETAP[K]
    graw = GTEMPO.integrate(LAT, GTEMPO._mult_A(GTEMPO.GTerm(BARP[K], b1; coeff=1), mpsK), mpsI)
    return z, 1 - real(graw / z)
end

function main()
    k = K
    npar = 2k*k + 1
    eds = [0.0, 0.35, 0.7, 1.05, 1.4]
    zanch = [fock_anchor(CELLS, RLINKS, e)[1] for e in eds]
    nanch = [real(fock_anchor(CELLS, RLINKS, e)[2]) for e in eds]
    @printf("锚 Z: %s\n", [round(v, digits=6) for v in zanch])
    @printf("锚 n: %s\n", [round(v, digits=6) for v in nanch])

    function resid(x; cache::Bool=false)
        out = ComplexF64[]
        for (i, e) in enumerate(eds)
            mpsK = build_kchain(e)
            z, n = composite(x, e)
            push!(out, ComplexF64(real(z) - zanch[i], 0.0))
            push!(out, ComplexF64(n - nanch[i], 0.0))
        end
        return out
    end

    Random.seed!(9)
    best = (Inf, nothing)
    for trial in 1:4
        x = vcat(0.2 .* randn(2k*k), [1.0])
        for it in 1:25
            F = resid(x)
            J = zeros(ComplexF64, length(F), npar)
            for p in 1:npar
                xp = copy(x); xp[p] += 1e-5
                J[:, p] = (resid(xp) - F) ./ 1e-5
            end
            dx = real.(pinv(J' * J + 1e-8 * I; rtol=1e-9) * (-J' * F))
            x .+= dx
            nf = norm(resid(x))
            it % 5 == 0 && (@printf("  trial %d iter %2d: |F| = %.3e\n", trial, it, nf); flush(stdout))
            nf < 1e-11 && break
            norm(dx) < 1e-12 && break
        end
        rf = norm(resid(x))
        rf < best[1] && (best = (rf, copy(x)))
        @printf("trial %d 完成: |F| = %.3e\n", trial, rf); flush(stdout)
        best[1] < 1e-11 && break
    end
    @printf("\n最优 |F| = %.3e\n", best[1])
    if best[1] < 1e-9
        x = best[2]
        println("端到端验证：")
        mpsI = build_imps(x)
        for e in eds
            mpsK = build_kchain(e)
            z = GTEMPO.integrate(LAT, mpsK, mpsI) * x[end]
            za, na = fock_anchor(CELLS, RLINKS, e)
            b1 = ETAP[K]
            graw = GTEMPO.integrate(LAT, GTEMPO._mult_A(GTEMPO.GTerm(BARP[K], b1; coeff=1), mpsK), mpsI)
            npipe = 1 - real(graw / z)
            @printf("εd=%.2f : Z err %.2e  n err %.2e\n", e, abs(real(z) - real(za)), abs(npipe - real(na)))
        end
    end
end
main()
