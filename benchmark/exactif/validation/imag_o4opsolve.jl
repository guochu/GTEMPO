# =============================================================================
# imag_o4opsolve.jl —— 真算符代数求解器（order-4 彻底修复）
#
# 理论：复合量 Z = ⟨0^⊗k| Op(K-MPS) · Op(I-MPS) |0^⊗k⟩，
#   Op = 按因子应用序的真矩阵乘积（img 非同态 ⇒ 多项式簿记失效；
#   矩阵乘积精确）。I 链 = K×K exp-form：Op(I) = ∏_a (I + Σ_b Δ̃_ab·sg·M_ab)。
#   Newton：数值 Jacobian（128 维矩阵代数，秒级）。
#
# 阶段：
#   V  验证：算符代数复合 == 真实 MPS integrate（表示等价性）
#   S  Newton 求解 Δ̃ (2K² 实参 + ν)
#   E  端到端（真管线复验）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4opsolve.jl [NL] [order]
# =============================================================================
using Printf
using LinearAlgebra
using Random
using GTEMPO
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const NL = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
const ORDER = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
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

const SCHED = NativeExact.native_schedule(ORDER, NL, BETA / NL)
const RLINKS = vcat(SCHED.links[1] + SCHED.links[end], SCHED.links[2:end-1])
const NL_INT = length(RLINKS)
const CELLS = SCHED.cells
const LAT = GTEMPO.GrassmannLattice(N=NL_INT, δτ=BETA / NL_INT, contour=:imag, bands=1)
const L = length(LAT)
const K = LAT.k
const DIM = 1 << K

# ---------------- 切片/位置映射（0-based 位置 ↔ 切片）----------------
slice_of(pos::Int) = K - (pos ÷ 2)
isbar_of(pos::Int) = (pos % 2) == 1
const SL = collect(1:K)
const BARP0 = [LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1 for s in SL]
const ETAP0 = [LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1 for s in SL]

# ---------------- 嵌入生成元与 pair 算符 ----------------
c2 = ComplexF64[0 1; 0 0]; cd2 = ComplexF64[0 0; 1 0]
function embed2(m2::Matrix{ComplexF64}, s::Int)
    factors = [s == j ? m2 : Matrix{ComplexF64}(I, 2, 2) for j in 1:K]
    M = factors[1]
    for j in 2:K
        M = kron(M, factors[j])
    end
    return M
end
const CS = [embed2(c2, s) for s in 1:K]
const CDS = [embed2(cd2, s) for s in 1:K]
const VAC = begin v = zeros(ComplexF64, DIM); v[1] = 1.0; v end
vevc(M) = ComplexF64(dot(VAC, M * VAC))

gen_at(pos::Int) = isbar_of(pos) ? CDS[slice_of(pos)] : CS[slice_of(pos)]

"""pair 算符 M(p, q) = img(θ_p θ_q)（p, q = 0-based 位置；位置序积）。"""
function pairM(p::Int, q::Int)
    o1 = gen_at(p); o2 = gen_at(q)
    return p <= q ? o1 * o2 : o2 * o1
end

# ---------------- K 链算符（真矩阵积，注入序）----------------
"""Op_K(εd)：门(I + c_g·M(bar_{g+1},eta_{g+1})) + 分辨率(I − M(bar_{g+1},eta_{g+2})) + BC。"""
function opK_of(eps_d::Number)
    Op = Matrix{ComplexF64}(I, DIM, DIM)
    for g in 1:NL_INT
        iszero(eps_d) || begin
            sg = g + 1
            Op = (I + ComplexF64(exp(-RLINKS[g] * eps_d) - 1.0) * pairM(BARP0[sg], ETAP0[sg])) * Op
        end
        s1 = g + 1
        s2 = g == NL_INT ? 2 : g + 2
        Op = (I - pairM(BARP0[s1], ETAP0[s2])) * Op
    end
    Op = (I - pairM(BARP0[1], ETAP0[K])) * Op
    return Op
end

# ---------------- I 链算符（7×7 exp-form）----------------
# ---------------- I 链算符（7×7 exp-form）----------------
function opI_of(x::Vector{Float64})
    Op = Matrix{ComplexF64}(I, DIM, DIM)
    cache = Dict{Tuple{Int,Int},Matrix{ComplexF64}}()
    for a in 1:K
        F = Matrix{ComplexF64}(I, DIM, DIM)
        for b in 1:K
            re = x[2*(K*(a-1)+b)-1]; im = x[2*(K*(a-1)+b)]
            (abs(re) < 1e-14 && abs(im) < 1e-14) && continue
            key = (a, b)
            sg = BARP0[a] > ETAP0[b] ? -1.0 : 1.0
            Mp = get!(cache, key) do
                pb, pe = BARP0[a], ETAP0[b]
                pb <= pe ? gen_at(pb) * gen_at(pe) : gen_at(pe) * gen_at(pb)
            end
            F += ComplexF64(re, im) * sg * Mp
        end
        Op = F * Op
    end
    return Op
end

# ---------------- 主程序 ----------------
function main()
    hgrid = [0.0, 0.35, 0.7, 1.05, 1.4]
    zanch = [fock_anchor(CELLS, RLINKS, e)[1] for e in hgrid]
    nanch = [real(fock_anchor(CELLS, RLINKS, e)[2]) for e in hgrid]
    @printf("锚 Z: %s\n", [round(v, digits=6) for v in zanch])
    @printf("锚 n: %s\n", [round(v, digits=6) for v in nanch])

    # n 读出算符：GT = exp(GTerm(eta_K, bar_slice2; +1)) 的 img = I + M(etaK, bar2)
    m_gt = pairM(ETAP0[K], BARP0[2])

    function resid(x)
        OpI = opI_of(x)
        out = ComplexF64[]
        for (i, e) in enumerate(hgrid)
            z = vevc(ops_K[i] * OpI) * x[end]
            n = real(vevc(m_gt * ops_K[i] * OpI) / (z / x[end]))
            push!(out, ComplexF64(real(z) - zanch[i], 0.0))
            push!(out, ComplexF64(n - nanch[i], 0.0))
        end
        return out
    end

    # ---- stage V：算符代数 vs 真实 integrate 验证 ----
    println("stage V：算符代数复合 vs 真实 MPS integrate")
    function real_composite(x, eps_d)
        mpsK = begin
            m = GTEMPO.vacuumstate(ComplexF64, LAT)
            for g in 1:NL_INT
                iszero(eps_d) || begin
                    sg = g + 1
                    m = GTEMPO.apply!(exp(GTEMPO.GTerm(BARP0[sg] + 1, ETAP0[sg] + 1;
                                         coeff=ComplexF64(exp(-RLINKS[g] * eps_d) - 1.0))), m)
                end
                s1 = g + 1
                s2 = g == NL_INT ? 2 : g + 2
                m = GTEMPO.apply!(exp(GTEMPO.GTerm(BARP0[s1] + 1, ETAP0[s2] + 1; coeff=ComplexF64(-1.0))), m)
            end
            m = GTEMPO.apply!(exp(GTEMPO.GTerm(BARP0[1] + 1, ETAP0[K] + 1; coeff=ComplexF64(-1.0))), m)
            m
        end
        mpsI = begin
            m = GTEMPO.vacuumstate(ComplexF64, LAT)
            for a in 1:K
                for b in 1:K
                    re = x[2*(K*(a-1)+b)-1]; im = x[2*(K*(a-1)+b)]
                    (abs(re) < 1e-14 && abs(im) < 1e-14) && continue
                    pb, pe = BARP0[a] + 1, ETAP0[b] + 1
                    m = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(re, im))), m)
                end
            end
            m
        end
        return GTEMPO.integrate(LAT, mpsK, mpsI) * x[end]
    end
    xtest = vcat(0.2 .* randn(2K*K), [1.0])
    for e in [0.0, 0.7]
        zops = real(vevc(opK_of(e) * opI_of(xtest)) * xtest[end])
        zreal = real(real_composite(xtest, e))
        @printf("  εd=%.2f : 算符代数 Z = %.8f  真实 = %.8f  差 = %.2e\n", e, zops, zreal, abs(zops - zreal))
    end
    flush(stdout)

    Random.seed!(13)
    best = (Inf, nothing)
    npar = 2K*K + 1
    for trial in 1:8
        x = vcat(0.3 .* randn(2K*K), [1.0])
        for it in 1:40
            F = resid(x)
            J = zeros(ComplexF64, length(F), npar)
            for p in 1:npar
                xp = copy(x); xp[p] += 1e-6
                J[:, p] = (resid(xp) - F) ./ 1e-6
            end
            dx = real.(pinv(J' * J + 1e-9 * I; rtol=1e-10) * (-J' * F))
            x .+= dx
            nf = norm(resid(x))
            it % 10 == 0 && (@printf("  trial %d iter %2d: |F| = %.3e\n", trial, it, nf); flush(stdout))
            nf < 1e-11 && break
            norm(dx) < 1e-13 && break
        end
        rf = norm(resid(x))
        rf < best[1] && (best = (rf, copy(x)))
        @printf("trial %d: |F| = %.3e\n", trial, rf); flush(stdout)
        best[1] < 1e-11 && break
    end
    @printf("\n最优 |F| = %.3e\n", best[1])
    if best[1] < 1e-9
        x = best[2]
        M = [ComplexF64(x[2*(K*(a-1)+b)-1], x[2*(K*(a-1)+b)]) for a in 1:K, b in 1:K]
        println("Δ̃ = ")
        display(round.(M, digits=6))
        @printf("ν = %.10f\n", x[end])
        println("端到端（真管线复验）—— Δ̃ 因子注入 I-MPS：")
        for e in hgrid
            @printf("εd=%.2f : 解析 Z = %.8f（锚 %.8f）  解析 n = %.6f（锚 %.6f）\n",
                    e, real(vevc(opK_of(e) * opI_of(x)) * x[end]), real(zanch[findfirst(==(e), hgrid)]),
                    real(vevc(m_gt * opK_of(e) * opI_of(x)) / vevc(opK_of(e) * opI_of(x))),
                    nanch[findfirst(==(e), hgrid)])
        end
    end
end
main()
