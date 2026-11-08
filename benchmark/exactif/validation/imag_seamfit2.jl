# =============================================================================
# imag_seamfit2.jl —— 干净行上的缝位 3×3 影响矩阵族拟合（闭式 Δ̃ 提取）
#
# 上一轮 seam_solve 败于污染行（g 混基 + z=0 守卫清零）。本轮：combo-3 K 链
# （c_gate = e^{−τε}−1, s_res = −1）+ h/g_raw 行（均 I 侧 Möbius 到单项式基）
# + 3×3 exp-form 族（Δ̃ 18 实参 + ν）多起点牛顿。命中 ⇒ Δ̃ = 闭式矩阵。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_seamfit2.jl
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const N = 2
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]

function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op_f(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = NativeExact._coupled_h(VB, HB)
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
annih_op_f(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end

const NL = 2
const L = 2 * NL + 2
const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
const EVENS = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const CFGIDX = [1 + sum(1 << p for p in A; init=0) for A in EVENS]
pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]
all_subsets(x) = isempty(x) ? [Tuple{Int,Int}[]] :
    [x[[m == 1 for m in msk]] for msk in Iterators.product(fill(0:1, length(x))...)]
idxof(A) = findfirst(==(A), EVENS)

function injected_mps(pairs::Vector{Tuple{Int,Int}})
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for (p1, p2) in sort(pairs)
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mps)
    end
    return mps
end

"""pipeline-v2 K 链（combo 3：c_gate = e^{−τε}−1，s_res = −1）。"""
function kchain_v2(eps_d)
    sched = NativeExact.native_schedule(2, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for g in 1:NL
        iszero(eps_d) && continue
        s = g + 1
        pb = LAT[GrassmannContourIndex(s, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe;
                             coeff=ComplexF64(exp(-ring_links[g] * eps_d) - 1.0))), mps)
    end
    for p in 1:N
        s1 = p + 1
        s2 = p == N ? 2 : p + 2
        pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(-1.0))), mps)
    end
    p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
    p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
    mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mps)
    return mps
end

# 3×3 切片位（切片 1=缝, 2, 3）
const SL = [1, 2, 3]
const BARP = [LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1 for s in SL]
const ETAP = [LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1 for s in SL]
const PR = [(BARP[a], ETAP[b], BARP[a] > ETAP[b] ? -1.0 : 1.0) for a in 1:3, b in 1:3]

"""3×3 exp-form IF：∏_a exp(Σ_b Δ̃[a,b]·θ_bar_a θ_eta_b)（行序 1,2,3）。"""
function if3_table(x)
    M = [ComplexF64(x[2*(3*(a-1)+b)-1], x[2*(3*(a-1)+b)]) for a in 1:3, b in 1:3]
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for a in 1:3
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        for b in 1:3
            iszero(M[a, b]) && continue
            pb, pe, sg = PR[a, b]
            push!(monos, (sort([pb, pe]), M[a, b] * sg))
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end

function main()
    hgrid = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875]
    ggrid = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 1.875]

    # ---- h 行（Möbius 到单项式基）
    hrows = Vector{Vector{ComplexF64}}()
    for ed in hgrid
        mpsK = kchain_v2(ed)
        gB = [GTEMPO.integrate(LAT, mpsK, injected_mps(pairs_of(B))) for B in EVENS]
        h = zeros(ComplexF64, length(EVENS))
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            for Bsub in all_subsets(pairs_of(B))
                sb = sort(collect(Iterators.flatten(Bsub)))
                sb = isempty(sb) ? Int[] : collect(sb)
                sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                val += sgn * gB[idxof(sb)]
            end
            h[ib] = val
        end
        push!(hrows, h)
        GC.gc(); flush(stdout)
    end
    # ---- g_raw 行（_mult_A 直测 + I 侧 Möbius）
    b1p = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
    grows = Vector{Vector{ComplexF64}}()
    for ed in ggrid
        mpsK = kchain_v2(ed)
        t_end = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1p; coeff=1)
        t_1 = GTEMPO.GTerm(LAT[GrassmannContourIndex(1, conj=false, branch=:τ, band=1)], b1p; coeff=1)
        mpsK_end = GTEMPO._mult_A(t_end, mpsK)
        mpsK_1 = GTEMPO._mult_A(t_1, mpsK)
        row = zeros(ComplexF64, 3 * length(EVENS))
        for (col, mpsKA) in ((0, mpsK), (length(EVENS), mpsK_end), (2length(EVENS), mpsK_1))
            raw = [GTEMPO.integrate(LAT, mpsKA, injected_mps(pairs_of(B))) for B in EVENS]
            for (ib, B) in enumerate(EVENS)
                val = 0.0 + 0.0im
                for Bsub in all_subsets(pairs_of(B))
                    sb = sort(collect(Iterators.flatten(Bsub)))
                    sb = isempty(sb) ? Int[] : collect(sb)
                    sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                    val += sgn * raw[idxof(sb)]
                end
                row[col + ib] = val
            end
        end
        push!(grows, row)
        GC.gc(); flush(stdout)
    end

    # ---- 联合行
    rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]
    for (r, ed) in enumerate(hgrid)
        push!(rows, hrows[r]); push!(rhs, fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)[1])
    end
    for (r, ed) in enumerate(ggrid)
        hz = hrows[findfirst(==(ed), hgrid)]
        vend = grows[r][length(EVENS) .+ (1:length(EVENS))]
        nanc = real(fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)[2])
        push!(rows, hz .- vend .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
    end
    A = hcat(rows...)
    @printf("联合行数 = %d\n", length(rhs))

    function resid(x; realΔ::Bool)
        itab = realΔ ? if3_table(vcat(x[1:9], zeros(9), x[10])) : if3_table(x)
        ivec = itab[CFGIDX]
        return [x[realΔ ? 10 : 19] * dot(rows[r], ivec) - rhs[r] for r in eachindex(rows)]
    end
    function newton(x0; iters=400, realΔ::Bool=false)
        x = copy(x0)
        npar = realΔ ? 10 : 19
        for it in 1:iters
            F = resid(x; realΔ=realΔ)
            J = zeros(ComplexF64, length(F), npar)
            for p in 1:npar
                xp = copy(x); xp[p] += 1e-7
                xm = copy(x); xm[p] -= 1e-7
                J[:, p] = (resid(xp; realΔ=realΔ) - resid(xm; realΔ=realΔ)) ./ 2e-7
            end
            dx = real.(pinv(J' * J + 1e-12 * I, rtol=1e-10) * (-J' * F))
            x += dx
            norm(F) < 1e-13 && break
            norm(dx) < 1e-16 && break
        end
        return x, norm(resid(x; realΔ=realΔ))
    end
    for realΔ in (true, false)
        Random.seed!(23)
        best = (Inf, nothing); hits = 0
        npar = realΔ ? 9 : 18
        for trial in 1:40
            x0 = vcat(0.6 .* randn(npar), [0.5 + rand()])
            x, rf = newton(x0; realΔ=realΔ)
            rf < best[1] && (best = (rf, copy(x)))
            rf < 1e-11 && (hits += 1)
            hits >= 3 && break
        end
        @printf("[realΔ=%d] 3×3 族命中 = %d，最优 |F| = %.3e\n", realΔ, hits, best[1])
        if best[1] < 1e-9
            global bestreal = best[2]
            x = best[2]
            M = [ComplexF64(x[2*(3*(a-1)+b)-1], x[2*(3*(a-1)+b)]) for a in 1:3, b in 1:3]
            println("\nΔ̃（行=切片 1缝,2,3 的 bar；列=切片 1缝,2,3 的 eta）= ")
            display(round.(M, digits=6))
            @printf("ν = %.10f\n", x[realΔ ? 10 : 19])
            realΔ && break
        end
    end
        println("\n对照：环 M 杂质块 / Schur（打印供解读）")
        h = NativeExact._coupled_h(VB, HB)
        mm = size(h, 1)
        cells = fill(BETA / N, N)
        Us = [exp(-cells[r] .* h) for r in 1:N]
        Mring = zeros(ComplexF64, mm * N, mm * N)
        for r in 1:N
            sgn = r == N ? +1.0 : -1.0
            for s2 in 1:mm, s1 in 1:mm
                Mring[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
            end
            nxt = r == N ? 1 : r + 1
            for s in 1:mm
                Mring[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0
            end
        end
        Iset = [(r-1)*mm + 1 for r in 1:N]
        Bset = [(r-1)*mm + s for r in 1:N for s in 2:mm]
        Dschur = -(Mring[Iset, Iset] - Mring[Iset, Bset] * (Mring[Bset, Bset] \ Mring[Bset, Iset]))
        println("环 M 杂质块（N×N）= ")
        display(round.(Mring[Iset, Iset], digits=6))
    println("−D_Schur = ")
    display(round.(-Dschur, digits=6))
end
main()
