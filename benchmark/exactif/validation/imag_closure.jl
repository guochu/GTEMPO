# =============================================================================
# imag_closure.jl —— pipeline-v2 闭环：I* 求解 + 结构解读 + 端到端 Z/n 验证
#
# combo 3（c_gate = e^{−τε}−1, s_res = −1）联合一致（残差 5e-11）⇒ 存在 I 多项式
# 同时满足 Z+n。本脚本：min-norm I* → Möbius pair-exp 实现 → 真实 GMPS 管线
# （注入 K 链 + 注入 I 链 + integrate + gf）端到端 vs Fock 锚（机器精度判据）。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_closure.jl [combo]
# =============================================================================
using Printf
using LinearAlgebra
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
pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]
all_subsets(x) = isempty(x) ? [Tuple{Int,Int}[]] :
    [x[[m == 1 for m in msk]] for msk in Iterators.product(fill(0:1, length(x))...)]
idxof(A) = findfirst(==(A), EVENS)

function injected_mps(pairs; coeffs=nothing)
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for (i, (p1, p2)) in enumerate(sort(pairs))
        c = coeffs === nothing ? 1.0 : coeffs[i]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(c))), mps)
    end
    return mps
end

function kchain_v2(eps_d; c_gate::Function, s_res::Float64)
    sched = NativeExact.native_schedule(2, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for g in 1:NL
        iszero(eps_d) && continue
        s = g + 1
        pb = LAT[GrassmannContourIndex(s, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe;
                             coeff=ComplexF64(c_gate(ring_links[g] * eps_d)))), mps)
    end
    for p in 1:N
        s1 = p + 1
        s2 = p == N ? 2 : p + 2
        pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(s_res))), mps)
    end
    p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
    p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
    mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mps)
    return mps
end

function main()
    combo = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3
    gate_c, s_res = [(τ -> exp(-τ) - 1.0, -1.0),
                     (τ -> 1.0 - exp(-τ), +1.0),
                     (τ -> exp(-τ) - 1.0, +1.0)][combo]
    println("combo $(combo)：c_gate = e^{-τε} − 1 系（s_res = $(s_res)）")

    hgrid = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875]
    ggrid = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 1.875]

    hrows = Vector{Vector{ComplexF64}}()
    for ed in hgrid
        mpsK = kchain_v2(ed; c_gate=gate_c, s_res=s_res)
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
    grows = Vector{Vector{ComplexF64}}()
    for ed in ggrid
        mpsK = kchain_v2(ed; c_gate=gate_c, s_res=s_res)
        row = zeros(ComplexF64, 3 * length(EVENS))
        b1p = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
        t_end = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1p; coeff=1)
        t_1 = GTEMPO.GTerm(LAT[GrassmannContourIndex(1, conj=false, branch=:τ, band=1)], b1p; coeff=1)
        mpsK_end = GTEMPO._mult_A(t_end, mpsK)
        mpsK_1 = GTEMPO._mult_A(t_1, mpsK)
        for (ib, B) in enumerate(EVENS)
            mpsI = injected_mps(pairs_of(B))
            row[ib] = hrows[findfirst(==(ed), hgrid)][ib]
            row[length(EVENS) + ib] = GTEMPO.integrate(LAT, mpsK_end, mpsI)
            row[2length(EVENS) + ib] = GTEMPO.integrate(LAT, mpsK_1, mpsI)
        end
        # I 侧 Möbius：pair-乘积注入值 → 纯单项式系数（与 h-行同基）
        gmono = zeros(ComplexF64, length(EVENS))
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            for ch in (:z, :ge, :g1)
                col = ch === :z ? 0 : ch === :ge ? length(EVENS) : 2length(EVENS)
                v = 0.0 + 0.0im
                for Bsub in all_subsets(pairs_of(B))
                    sb = sort(collect(Iterators.flatten(Bsub)))
                    sb = isempty(sb) ? Int[] : collect(sb)
                    sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                    v += sgn * row[col + idxof(sb)]
                end
                # 写回时只对当前通道操作——改为分别处理
            end
        end
        # （逐通道 Möbius）
        for (col, off) in ((0, 0), (length(EVENS), length(EVENS)), (2length(EVENS), 2length(EVENS)))
            for (ib, B) in enumerate(EVENS)
                v = 0.0 + 0.0im
                for Bsub in all_subsets(pairs_of(B))
                    sb = sort(collect(Iterators.flatten(Bsub)))
                    sb = isempty(sb) ? Int[] : collect(sb)
                    sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                    v += sgn * row[col + idxof(sb)]
                end
                row[off + ib] = v
            end
        end
        push!(grows, row)
        GC.gc(); flush(stdout)
    end

    # 联合行（Z + n，1−g_end 变体）
    rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]
    for (r, ed) in enumerate(hgrid)
        push!(rows, hrows[r]); push!(rhs, fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)[1])
    end
    for (r, ed) in enumerate(ggrid)
        hz = hrows[findfirst(==(ed), hgrid)]
        ge = grows[r][length(EVENS) .+ (1:length(EVENS))]
        nanc = real(fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)[2])
        push!(rows, hz .- ge .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
    end
    A = hcat(rows...)
    Ist = A * ((A' * A + 1e-13 * I) \ rhs)
    @printf("联合残差 = %.3e\n", norm(A' * Ist - rhs))
    println("\nI* 非零项：")
    for (i, Aset) in enumerate(EVENS)
        v = Ist[i]
        abs(v) > 1e-8 || continue
        @printf("   %-16s : %11.7f%s\n", Aset, real(v),
                abs(imag(v)) > 1e-8 ? " + $(round(imag(v), digits=7))im" : "")
    end

    # 全偶单项式序贯消元实现：I = ∏_S exp(c_S·θ_S)（GTerm 支持任意偶腿）
    println("\n── 全单项式序贯消元实现 ──")
    ν_eff = Ist[1]                       # min-norm 解的整体标度
    println("ν_eff（I* ∅ 系数）= ", ν_eff)
    ivec = zeros(ComplexF64, 1 << L)
    for (i, Aset) in enumerate(EVENS)
        ivec[1 + sum(1 << p for p in Aset; init=0)] = Ist[i] / ν_eff
    end
    factors = Vector{Tuple{Vector{Int},ComplexF64}}()
    tcur = zeros(ComplexF64, 1 << L); tcur[1] = 1.0
    for Aset in sort(EVENS; by=length)
        want = ivec[1 + sum(1 << p for p in Aset; init=0)]
        isempty(Aset) && (abs(want - 1) > 1e-9 && error("∅ 系数 ≠ 1"))
        isempty(Aset) && continue
        c = want - tcur[1 + sum(1 << p for p in Aset; init=0)]
        abs(c) < 1e-10 && continue
        push!(factors, (Aset, c))
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im), (Aset, c)]
        tcur = mul_monomials!(zeros(ComplexF64, 1 << L), tcur, L, Monomials(monos))
    end
    err = maximum(abs.(tcur - ivec))
    @printf("实现因子数 = %d，重构误差 = %.3e\n", length(factors), err)
    err < 1e-9 || error("序贯消元失败")
    for (Aset, c) in factors
        @printf("   %-16s : %.7f%s\n", Aset, real(c),
                abs(imag(c)) > 1e-9 ? " + $(round(imag(c), digits=7))im" : "")
    end

    # 端到端：注入 I*（逐因子 exp(GTerm)）+ 各 εd 的 Z/n vs 锚
    println("\n── 端到端（真实 GMPS 注入链 + integrate + gf）──")
    mpsI_build = () -> begin
        m = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (Aset, c) in factors
            m = GTEMPO.apply!(exp(GTEMPO.GTerm([p + 1 for p in Aset]...; coeff=ComplexF64(c))), m)
        end
        return m
    end
    # —— 诊断：复合 mpsI 的多项式 vs ivec（单腿提取器）
    mps_to_dense = (x) -> begin
        A = ones(ComplexF64, 1, 1)
        for j in 1:length(x)
            t = x[j]
            spL, domR = codomain(t), domain(t)
            dL, dR = spL[1], domR[1]
            dimL = [blockdim(dL, Z2Irrep(q)) for q in 0:1]
            dimR = [blockdim(dR, Z2Irrep(q)) for q in 0:1]
            newA = zeros(ComplexF64, dimR[1] + dimR[2], 2 * size(A, 2))
            for (f1, f2) in fusiontrees(t)
                qL = f1.uncoupled[1].n; bp = f1.uncoupled[2].n; qR = f2.uncoupled[1].n
                mm = t[f1, f2]
                for α in 1:dimL[qL+1], β in 1:dimR[qR+1], cfg in 1:size(A, 2)
                    newA[(qR == 0 ? β : dimR[1] + β), cfg + (bp << (j - 1))] +=
                         mm[α, 1, β] * A[(qL == 0 ? α : dimL[1] + α), cfg]
                end
            end
            A = newA
        end
        vec((GTEMPO.scaling(x)^L) .* A)
    end
    mI = mpsI_build()
    ivec2 = zeros(ComplexF64, 1 << L)
    for (i, Aset) in enumerate(EVENS)
        ivec2[1 + sum(1 << p for p in Aset; init=0)] = Ist[i] / ν_eff
    end
    extI = mps_to_dense(mI)
    dmax = maximum(abs.(extI - ivec2))
    println("诊断：复合 mpsI 提取 vs I*/ν_eff max|Δ| = ", dmax)
    if dmax > 1e-10
        for idx in 1:length(ivec2)
            abs(extI[idx] - ivec2[idx]) > 1e-10 || continue
            bits = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
            @printf("   %-16s : 提取 %10.6f vs 真值 %10.6f\n", bits, real(extI[idx]), real(ivec2[idx]))
        end
    end
    for ed in [0.0, 0.35, 0.7, 1.05, 1.4]
        mpsK = kchain_v2(ed; c_gate=gate_c, s_res=s_res)
        mpsI = mpsI_build()
        z0 = GTEMPO.integrate(LAT, mpsK, mpsI)
        z = z0 * ν_eff
        # —— gf 线性判决：g_raw(复合) vs Σ I_b·g_raw(P_b)
        b1p = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
        tend = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1p; coeff=1)
        gcomp = GTEMPO.integrate(LAT, GTEMPO._mult_A(tend, mpsK), mpsI)
        gcomp_raw = gcomp
        # Σ I_b·g_raw(P_b)（用 Ist 系数；g_raw(P_b) 重测）
        gsum = 0.0 + 0.0im
        for (i, Aset) in enumerate(EVENS)
            abs(Ist[i]) > 1e-12 || continue
            mpsIb = injected_mps(pairs_of(Aset))
            gsum += Ist[i] * GTEMPO.integrate(LAT, GTEMPO._mult_A(tend, mpsK), mpsIb)
        end
        @printf("  [线性判决 εd=%.2f] g_raw(复合) = %10.6f  Σ I·g_raw = %10.6f  差 = %9.2e\n",
                ed, real(gcomp_raw), real(gsum), abs(gcomp_raw - gsum))
        gend = gcomp / z0
        za, na = fock_anchor(fill(BETA / N, N), fill(BETA / N, N), ed)
        @printf("εd=%.2f : Z = %12.8f（锚 %12.8f, err %9.2e）  1−g_end = %9.6f（n锚 %9.6f, err %9.2e）\n",
                ed, real(z), real(za), abs(real(z) - real(za)), 1 - real(gend), real(na), abs(1 - real(gend) - real(na)))
    end
end
main()
