# =============================================================================
# imag_lintest.jl —— canonical 线性判决：稀疏解的探针多项式叠加 vs canonical 预测
#
# 问题：复合量 canonical 线性（composite(K,I) = Σ_b I_b·h_b(K)）对构建型 MPS
#   是否成立？若成立 → 单 MPS（序贯消元构建）路线可行；否则 → 探针循环为最终形态。
# 方法：OMP 稀疏解 → 每探针提 canonical 多项式（mps_to_dense）→ I_b = Σ c_j·S_b
#   → canonical 预测 Σ_b I_b·h_b(εd) vs fock_anchor 真值。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_lintest.jl
# =============================================================================
using Printf
using Random
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation, Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const WB = parse(Float64, get(ENV, "IMAG_WB", "1.0"))
const AB = parse(Float64, get(ENV, "IMAG_AB", "0.5"))
const TAG = get(ENV, "IMAG_TAG", "n2o4ma")
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

function mps_to_dense(x::GTEMPO.GrassmannMPS, L::Int)
    n = length(x)
    n == L || throw(ArgumentError("site 数应等于 L"))
    A = ones(ComplexF64, 1, 1)
    for j in 1:n
        t = x[j]
        spL, domR = codomain(t), domain(t)
        dL, dR = spL[1], domR[1]
        dimL = [blockdim(dL, Z2Irrep(q)) for q in 0:1]
        dimR = [blockdim(dR, Z2Irrep(q)) for q in 0:1]
        offR = (0, dimR[1])
        newA = zeros(ComplexF64, dimR[1] + dimR[2], 2 * size(A, 2))
        for (f1, f2) in fusiontrees(t)
            qL = f1.uncoupled[1].n; bp = f1.uncoupled[2].n; qR = f2.uncoupled[1].n
            m = t[f1, f2]
            for α in 1:dimL[qL+1], β in 1:dimR[qR+1], cfg in 1:size(A, 2)
                newA[offR[qR+1] + β, cfg + (bp << (j - 1))] +=
                     m[α, 1, β] * A[α + (qL == 1 ? dimL[1] : 0), cfg]
            end
        end
        A = newA
    end
    return vec((scaling(x)^L) .* A)
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    NL, order = 2, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    cells = sched.cells

    # 稀疏解
    sol = Vector{Tuple{String,Int,ComplexF64}}()
    for (i, ln) in enumerate(eachline(joinpath(outdir, "$(TAG)_dsol_omp.csv")))
        i == 1 && continue
        f = split(strip(ln), ",")
        push!(sol, (f[1], parse(Int, f[2]), ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))))
    end
    println("稀疏解项数：", length(sol)); flush(stdout)

    inj = function(A::Vector{Int}, σ::Vector{Float64})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (i, c) in zip(1:2:length(A), σ)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(A[i] + 1, A[i+1] + 1; coeff=ComplexF64(c))), mps)
        end
        return mps
    end

    # 探针 canonical 多项式
    rng = MersenneTwister(20261109)
    sprobes = Vector{Tuple{Vector{Int}, Vector{Float64}}}()
    for B in EVENS
        np2 = length(B) ÷ 2
        push!(sprobes, (B, [rand(rng) < 0.5 ? -1.0 : 1.0 for _ in 1:np2]))
    end
    pprobes = Vector{Tuple{Vector{Int}, Vector{Int}}}()
    for _ in 1:NE
        push!(pprobes, (EVENS[rand(rng, 1:NE)], EVENS[rand(rng, 1:NE)]))
    end

    Ipoly = zeros(ComplexF64, 1 << L)
    for (k, (fam, jp, c)) in enumerate(sol)
        mps = fam == "inj" ? inj(EVENS[jp], ones(Float64, length(EVENS[jp]) ÷ 2)) :
              fam == "s" ? inj(sprobes[jp][1], sprobes[jp][2]) : begin
                  (B1, B2) = pprobes[jp]
                  m1 = inj(B1, ones(Float64, length(B1) ÷ 2))
                  m2 = inj(B2, ones(Float64, length(B2) ÷ 2))
                  GTEMPO.mult!(m1, m2; trunc=NoTruncation())[1]
              end
        Sp = mps_to_dense(mps, L)
        Ipoly .+= c .* Sp
        k % 50 == 0 && (println("  #$k/$(length(sol))"); flush(stdout))
    end
    println("I_poly 支集：", count(x -> abs(x) > 1e-12, Ipoly)); flush(stdout)

    # canonical 预测：Z_can = Σ_b I_b·h_b(εd)
    hrows = Dict{Float64, Vector{ComplexF64}}()
    for file in sort(filter(f -> occursin("imag_$(TAG)_h", f) && endswith(f, ".csv"), readdir(outdir, join=true)))
        for (i, ln) in enumerate(eachline(file))
            i == 1 && continue
            f = split(strip(ln), ",")
            get!(hrows, parse(Float64, f[1]), zeros(ComplexF64, NE))[parse(Int, f[2])] =
                ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
        end
    end
    println("\ncanonical 线性判决（Z_can = Σ_b I_b·h_b vs 锚）：")
    for ed in [0.0, 0.35, 0.7, 1.05, 1.4]
        h = hrows[ed]
        # h_b ↔ EVENS[ib]：Z_can = Σ_ib I_b(ib)·h[ib]
        zcan = 0.0 + 0.0im
        for (ib, B) in enumerate(EVENS)
            m = sum(1 << p for p in B; init=0)
            v = Ipoly[m + 1]
            abs(v) > 1e-14 || continue
            zcan += v * h[ib]
        end
        (za, na) = fock_anchor(cells, ring_links, ed)
        @printf("  εd=%-5.2f : Z_can = %12.6f（真值 %12.6f, 违反 %.2e）\\n",
                ed, real(zcan), real(za), abs(real(zcan) - real(za)))
        flush(stdout)
    end
end
main()
