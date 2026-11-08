# =============================================================================
# imag_identify.jl —— I* 的闭式识别：环杂质扇区候选 vs 数值解逐系数对照
#
# 候选（均映射 ring pos p → 管线切片 p+1）：
#   c1  T_imp(εd=0)：全 kernel+resolution 环多项式 bath 积掉
#   c2  c1 去杂质 resolution pair（K 链已带）
#   c3  exp-form(−D_Schur)（稀疏 exp-form，canonical）
# 目标：combo-3 K 链联合行（干净基）的 min-norm I*。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_identify.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations, permutations
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_dense.jl"))
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

"""imag_dense.dense_ring 的表版本（逐字移植）。"""
function dense_ring_table(Ws)
    Nn = length(Ws)
    nm = size(Ws[1], 1) |> x -> round(Int, log2(x))
    Ld = 2nm * Nn
    Ld > 22 && error("稠密规模过大")
    dense = zeros(ComplexF64, 1 << Ld)
    dense[1] = 1.0
    bar(p, s) = 2 * ((p - 1) * nm + s - 1)
    eta(p, s) = bar(p, s) + 1
    for p in 1:Nn
        kmonos = kernel_monomials([bar(p, s) for s in 1:nm],
                                  [eta(p, s) for s in 1:nm], Ws[p])
        if p == Nn
            kmonos = [(loc, isodd(count(b -> b in [bar(Nn, s2) for s2 in 1:nm], loc)) ?
                       -c : c) for (loc, c) in kmonos]
        end
        dense = mul_monomials!(zeros(ComplexF64, 1 << Ld), dense, Ld, Monomials(kmonos))
        nxt = p == Nn ? 1 : p + 1
        for s in 1:nm
            rmonos = Monomials([(Int[], 1.0 + 0.0im),
                                ([bar(p, s), eta(nxt, s)], -1.0 + 0.0im)])
            dense = mul_monomials!(zeros(ComplexF64, 1 << Ld), dense, Ld, rmonos)
        end
    end
    return dense, Ld
end

"""环 M 矩阵（全模）与杂质 Iset/Bset（ring_delta 内部结构）。"""
function ring_M()
    h = NativeExact._coupled_h(VB, HB)
    mm = size(h, 1)
    cells = fill(BETA / N, N)
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0
        end
    end
    Iset = [(r-1)*mm + 1 for r in 1:N]
    Bset = [(r-1)*mm + s for r in 1:N for s in 2:mm]
    return M, Iset, Bset, mm
end

function main()
    hgrid = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875]
    ggrid = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 1.875]
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
    Ist = A * ((A' * A + 1e-13 * I) \ rhs)
    @printf("联合残差 = %.3e\n", norm(A' * Ist - rhs))
    ivec = zeros(ComplexF64, 1 << L)
    for (i, Aset) in enumerate(EVENS)
        ivec[1 + sum(1 << p for p in Aset; init=0)] = Ist[i]
    end
    ν_eff = Ist[1]
    target = ivec ./ ν_eff
    println("\nI*/ν_eff 非零项：")
    for (i, Aset) in enumerate(EVENS)
        v = Ist[i] / ν_eff
        abs(v) > 1e-8 || continue
        @printf("   %-16s : %11.7f\n", Aset, real(v))
    end

    # ---- 候选构造（ring → 管线位映射；ring pos p 的 bar/eta ↔ 切片 p+1）
    rpos(p, isbar) = begin
        s = p + 1
        isbar ? LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1 :
                LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1
    end
    Ws, _ = interval_operators(fill(BETA / N, N), fill(BETA / N, N), 0.0)
    tot, Ld = dense_ring_table(Ws)
    bath_gens = Int[]
    for p in 1:N
        push!(bath_gens, 4 * (p - 1) + 2, 4 * (p - 1) + 3)
    end
    Timp = integrate_out(tot, Ld, bath_gens)
    ν_r = Timp[1]
    Tm = zeros(ComplexF64, 1 << L)
    for idx in 1:length(Timp)
        Timp[idx] == 0 && continue
        bits = [p for p in 0:(Ld-1) if (idx - 1) >> p & 1 == 1]
        nb = [rpos(b ÷ 4 + 1, b % 4 == 0) for b in bits]
        Tm[1 + sum(1 << b for b in nb; init=0)] += Timp[idx]
    end
    c1 = Tm ./ Tm[1]

    # c2：去杂质 resolution pair（pair(bar_p, eta_{p+1})，环闭合）——逐 pair exp 除
    Mcand = copy(c1)
    for p in 1:N
        pb = rpos(p, true); pe = rpos(p == N ? 1 : p + 1, false)
        # 该 pair 的 canonical monomial 系数除以 −1（resolution factor 1 − θθ）
        key = 1 + (1 << pb) + (1 << pe)
        Mcand[key] = Mcand[key] + 1.0        # 去掉 −1 分辨率：+1 抵消
    end
    c2 = Mcand ./ Mcand[1]

    # c3：exp-form(−D_Schur)
    M, Iset, Bset, mm = ring_M()
    Dschur = -(M[Iset, Iset] - M[Iset, Bset] * (M[Bset, Bset] \ M[Bset, Iset]))
    c3 = zeros(ComplexF64, 1 << L); c3[1] = 1.0
    for j in 1:N, i in 1:N
        iszero(Dschur[i, j]) && continue
        pb = rpos(i, true); pe = rpos(j, false)
        sg = pb > pe ? -1.0 : 1.0
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im),
                                               (sort([pb, pe]), -Dschur[i, j] * sg)]
        c3 = mul_monomials!(zeros(ComplexF64, 1 << L), c3, L, Monomials(monos))
    end
    c3 = c3 ./ c3[1]

    # ---- F2：有序实行积族拟合（直接对 target，无 oracle）
    # 行 a = ∏_b exp(Δ_ab·sg_ab·θθ)，Δ 实 3×3，行序 π；参数 = Δ(9 实) + ln C
    SL = [1, 2, 3]
    BARP3 = [LAT[GrassmannContourIndex(sl, :τ, 1, true)] - 1 for sl in SL]
    ETAP3 = [LAT[GrassmannContourIndex(sl, :τ, 1, false)] - 1 for sl in SL]
    PR3 = [(BARP3[a], ETAP3[b], BARP3[a] > ETAP3[b] ? -1.0 : 1.0) for a in 1:3, b in 1:3]
    function fam_table(x, order)
        M = reshape(x[1:9], 3, 3)
        acc = zeros(ComplexF64, 1 << L); acc[1] = exp(x[10])
        for a in order
            for b in 1:3
                iszero(M[a, b]) && continue
                pb, pe, sg = PR3[a, b]
                monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im),
                                                       (sort([pb, pe]), M[a, b] * sg)]
                acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
            end
        end
        return acc
    end
    bestfam = (Inf, nothing, nothing)
    for perm in permutations((1, 2, 3))
        # 多起点牛顿（实 10 参数）
        for trial in 1:6
            x0 = vcat(0.5randn(9), [0.0])
            xx = copy(x0)
            for it in 1:200
                Fv = [fam_table(xx, perm)[CFGIDX[i]] - target[CFGIDX[i]] for i in 1:length(EVENS)]
                J = zeros(ComplexF64, length(Fv), 10)
                for p in 1:10
                    xp = copy(xx); xp[p] += 1e-6
                    J[:, p] = [fam_table(xp, perm)[CFGIDX[i]] - fam_table(xx, perm)[CFGIDX[i]] for i in 1:length(EVENS)] ./ 1e-6
                end
                dx = real.(pinv(J' * J + 1e-12 * I) * (-J' * Fv))
                xx += dx
                norm(Fv) < 1e-13 && break
            end
            rf = norm([fam_table(xx, perm)[CFGIDX[i]] - target[CFGIDX[i]] for i in 1:length(EVENS)])
            if rf < bestfam[1]
                bestfam = (rf, copy(xx), perm)
            end
            rf < 1e-12 && break
        end
        bestfam[1] < 1e-12 && break
    end
    @printf("\nF2 有序实行积：最优 |F| = %.3e（行序 %s）\n", bestfam[1], string(bestfam[3]))
    if bestfam[1] < 1e-10
        M = reshape(bestfam[2][1:9], 3, 3)
        println("Δ̃（实，行=切片1缝,2,3；列=eta1缝,2,3）= ")
        display(round.(M, digits=6))
        @printf("C = %.6f\n", exp(bestfam[2][10]))
    end

    println("\n候选 vs I*/ν_eff（|Δ| 最大）与 ∅：")
    for (nm, cand) in [("c1 T_imp 映射", c1), ("c2 c1 去分辨率", c2), ("c3 exp(−D_Schur)", c3)]
        err = maximum(abs.(cand - target))
        @printf("   %-18s : max err = %9.3e   ∅ = %.6f\n", nm, err, real(cand[1]))
    end
    println("\n目标（I*/ν_eff）pair 项 vs c1 pair 项：")
    for (i, Aset) in enumerate(EVENS)
        (length(Aset) == 2) || continue
        @printf("   %-8s : 目标 %10.6f  c1 %10.6f  c3 %10.6f\n", Aset,
                real(target[CFGIDX[i]]), real(c1[CFGIDX[i]]), real(c3[CFGIDX[i]]))
    end
end
main()
