# =============================================================================
# imag_trace_oracle.jl —— trace 泛函黑盒全测 + τ-稠密模拟 + 门情形验证
#
# 原理：integrate() 对 (K 链多项式, I 链多项式) 是双线性泛函
#       Z(S, I) = Σ_C (S⊗I)_canon[C] · τ[C]
# τ = 单项式对→值。偶集合的多项式 = ∏_pairs exp(θθ)（canonical +1），故
#   g(A,B) = Σ_{A'⊆pairs(A), B'⊆pairs(B)} τ[set(A')∪set(B')]
# Möbius 反演读全 τ（even set → 相邻连对分解唯一）。τ 完整后一切判决 =
# 稠密多项式代数（机器速度）。
#
# 阶段：
#   A  读 τ（even-even 块，32×32 = 1024 oracle 调用）
#   B1 探针组复验（τ-模拟 vs 真实 integrate）
#   B2 真骨架(εd=0)×真IF(Δ=D_full) 复验
#   B3 εd=0.7 门情形（解析 K 链 gate_sign ±标定 vs 真管线）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_trace_oracle.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))   # Monomials/mul_monomials!
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

# ------------------------------------------------------------------ 模型 --
const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
const EPSD = 0.7
hcoupled() = NativeExact._coupled_h(VB, HB)

function ring_delta(cells::Vector{Float64})
    N = length(cells)
    h = hcoupled()
    mm = size(h, 1)
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
    D = M[Iset, Iset] - M[Iset, Bset] * (M[Bset, Bset] \ M[Bset, Iset])
    Mb = zeros(ComplexF64, (mm-1)*N, (mm-1)*N)
    hbb = exp(-(BETA / N) .* HB)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm-1, s1 in 1:mm-1
            Mb[(r-1)*(mm-1) + s2, (r-1)*(mm-1) + s1] += sgn * hbb[s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm-1
            Mb[(r-1)*(mm-1) + s, (nxt-1)*(mm-1) + s] += 1.0
        end
    end
    return -D, det(M) / det(Mb)
end

mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
end
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex) = x.Δ[i.j, j.j]

"""vacuum + 注入 pair 集（canonical +1）。"""
function injected_mps(lat, pairs::Vector{Tuple{Int,Int}})
    mps = GTEMPO.vacuumstate(ComplexF64, lat)
    for (p1, p2) in sort(pairs)
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mps)
    end
    return mps
end

"""真实管线 Z（gates + BC + hybriddynamics + integrate·ν）。"""
function run_pipeline_real(Δ, ν, eps_d; N=2, order=2)
    sched = NativeExact.native_schedule(order, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(Δ); band=1, trunc=NoTruncation())
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    model = ToulouseIM(ϵ_d=eps_d)
    for g in 1:nl
        fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
        sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
        mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    return GTEMPO.integrate(lat, mpsK, mpsI) * ν
end

"""GMPS → canonical 稠密表（poly = scaling^L × 链）。"""
function mps_to_dense(x, L)
    n = length(x)
    n == L || error("n=L 预期（n=$n L=$L）")
    A = ones(ComplexF64, 1, 1)
    for j in 1:n
        t = x[j]
        spL, domR = codomain(t), domain(t)
        (length(spL) == 2 && length(domR) == 1) || error("site $j 结构异常")
        dL, dR = spL[1], domR[1]
        dimL = [blockdim(dL, Z2Irrep(q)) for q in 0:1]
        dimR = [blockdim(dR, Z2Irrep(q)) for q in 0:1]
        newA = zeros(ComplexF64, dimR[1] + dimR[2], 2 * size(A, 2))
        for (f1, f2) in fusiontrees(t)
            qL = f1.uncoupled[1].n; bp = f1.uncoupled[2].n; qR = f2.uncoupled[1].n
            m = t[f1, f2]
            for α in 1:dimL[qL+1], β in 1:dimR[qR+1], cfg in 1:size(A, 2)
                newA[(qR == 0 ? β : dimR[1] + β), cfg + (bp << (j - 1))] +=
                     m[α, 1, β] * A[(qL == 0 ? α : dimL[1] + α), cfg]
            end
        end
        A = newA
    end
    return vec((scaling(x)^L) .* A)
end

"""分析 IF 表：∏ᵢ exp(ā_i Σⱼ Δ_ij a_j)（canonical；与 hybriddynamics 已验证一致）。"""
function analytic_if(Δ, lat, L)
    Nτ = lat.Nτ
    barpos(i) = lat[GrassmannContourIndex(i + 1, :τ, 1, true)]
    etapos(j) = lat[GrassmannContourIndex(j + 1, :τ, 1, false)]
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for i in 1:Nτ
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        for j in 1:Nτ
            iszero(Δ[i, j]) && continue
            pb, pe = barpos(i) - 1, etapos(j) - 1
            sgn = (pb > pe) ? -1.0 : 1.0
            push!(monos, (sort([pb, pe]), Δ[i, j] * sgn))
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end

"""解析 K 链表（门 g@切片 g+1：kernel = exp(ξ̄(e^{−τεn̂}−I)ξ)，canonical
pair 系数 gate_sign·(e^{−τε_g}−1)；BC 最后乘）。εd=0 时应等于提取骨架。"""
function analytic_kchain(eps_d, nl, lat, L; gate_sign::Real=+1.0)
    sched = NativeExact.native_schedule(2, nl, BETA / nl)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for g in 1:nl
        iszero(eps_d) && continue
        s = g + 1
        pb = lat[GrassmannContourIndex(s, :τ, 1, true)] - 1
        pe = lat[GrassmannContourIndex(s, :τ, 1, false)] - 1
        c = gate_sign * (exp(-ring_links[g] * eps_d) - 1.0)
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L,
                             Monomials([(Int[], 1.0 + 0.0im),
                                        (sort([pb, pe]), ComplexF64(c))]))
    end
    pbc1 = lat[GrassmannContourIndex(1, :τ, 1, true)] - 1
    pbc2 = lat[GrassmannContourIndex(lat.k, :τ, 1, false)] - 1
    acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L,
                         Monomials([(Int[], 1.0 + 0.0im),
                                    (sort([pbc1, pbc2]), ComplexF64(1.0))]))
    return acc
end

"""canonical 位集多项式乘（表 → 表）。"""
function mul_tables(a::Vector{ComplexF64}, b::Vector{ComplexF64}, L::Int)
    out = zeros(ComplexF64, 1 << L)
    monos = Tuple{Vector{Int}, ComplexF64}[]
    for idx in 1:length(b)
        b[idx] == 0 && continue
        bits = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
        push!(monos, (bits, b[idx]))
    end
    return mul_monomials!(out, a, L, Monomials(monos))
end

# ================================================================== 主程序 --
function main()
    N = 2
    nl = 2
    cells = fill(BETA / N, N)
    sched = NativeExact.native_schedule(2, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(lat)

    # ---------------- A：读 τ ----------------
    evensets = Vector{Vector{Int}}([Int[]])
    for k in (2, 4, 6)
        for comb in combinations(0:(L-1), k)
            push!(evensets, collect(comb))
        end
    end
    ns = length(evensets)
    @printf("偶集合数 = %d（每链）\n", ns)

    pairs_of(A) = Tuple{Int,Int}[(A[2i-1], A[2i]) for i in 1:(length(A) ÷ 2)]
    all_subsets(x) = isempty(x) ? [Tuple{Int,Int}[]] :
        [x[[m == 1 for m in msk]] for msk in Iterators.product(fill(0:1, length(x))...)]

    println("oracle 网格 32×32 …")
    ggrid = Array{ComplexF64}(undef, ns, ns)
    for (ia, A) in enumerate(evensets), (ib, B) in enumerate(evensets)
        ggrid[ia, ib] = GTEMPO.integrate(lat, injected_mps(lat, pairs_of(A)),
                                         injected_mps(lat, pairs_of(B)))
        (ia * ns + ib) % 64 == 0 && (GC.gc(); flush(stdout))
    end
    GC.gc(); flush(stdout)

    idxof(A) = findfirst(==(A), evensets)
    # τ2[(A,B)]：K 侧单项式 A × I 侧单项式 B → 值（键 = 对，非并集！）
    τ2 = Dict{Tuple{Vector{Int},Vector{Int}}, ComplexF64}()
    for (ia, A) in enumerate(evensets), (ib, B) in enumerate(evensets)
        val = 0.0 + 0.0im
        for Asub in all_subsets(pairs_of(A)), Bsub in all_subsets(pairs_of(B))
            sa = sort(collect(Iterators.flatten(Asub))); sb = sort(collect(Iterators.flatten(Bsub)))
            sa = isempty(sa) ? Int[] : collect(sa); sb = isempty(sb) ? Int[] : collect(sb)
            sgn = (-1.0)^((length(A) - length(sa)) ÷ 2 + (length(B) - length(sb)) ÷ 2)
            val += sgn * ggrid[idxof(sa), idxof(sb)]
        end
        τ2[(A, B)] = val
    end
    nz = count(v -> abs(v) > 1e-12, values(τ2))
    @printf("τ2 条目 = %d（非零 %d）；τ2[(∅,∅)] = %s\n", length(τ2), nz, τ2[(Int[], Int[])])
    GC.gc(); flush(stdout)

    # ---------------- τ-稠密模拟 ----------------
    """Z_emul：S、I 表（L 位 canonical 系数）→ 直接双和 τ2（无排列符号）。"""
    function Zemul(Stab, Itab)
        configs(t) = begin
            out = Tuple{Vector{Int}, ComplexF64}[]
            for idx in 1:length(t)
                t[idx] == 0 && continue
                bits = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
                push!(out, (bits, t[idx]))
            end
            out
        end
        s = 0.0 + 0.0im
        for (a, ca) in configs(Stab), (b, cb) in configs(Itab)
            v = get(τ2, (a, b), 0.0)
            v == 0 && continue
            s += ca * cb * v
        end
        return s
    end

    # ---------------- B1：探针复验 ----------------
    println("\n── B1：探针组复验 ──")
    probes = [
        (K=Tuple{Int,Int}[], I=Tuple{Int,Int}[]),
        (K=[(0, 1)], I=Tuple{Int,Int}[]),
        (K=[(0, 1)], I=[(0, 1)]),
        (K=[(2, 3)], I=[(0, 1)]),
        (K=[(0, 3)], I=Tuple{Int,Int}[]),
        (K=[(0, 1), (2, 3)], I=[(4, 5)]),
        (K=[(0, 1)], I=[(0, 1), (2, 3)]),
        (K=[(0, 2)], I=[(1, 3)]),
        (K=[(1, 4)], I=Tuple{Int,Int}[]),
    ]
    for pr in probes
        mk = injected_mps(lat, collect(Tuple{Int,Int}, pr.K))
        mi = injected_mps(lat, collect(Tuple{Int,Int}, pr.I))
        zr = GTEMPO.integrate(lat, mk, mi)
        stab = zeros(ComplexF64, 1 << L); stab[1] = 1.0
        for (p1, p2) in sort(collect(Tuple{Int,Int}, pr.K))
            v = zeros(ComplexF64, 1 << L); v[1] = 1.0; v[1 + (1 << p1) + (1 << p2)] = 1.0
            stab = mul_tables(stab, v, L)
        end
        itab = zeros(ComplexF64, 1 << L); itab[1] = 1.0
        for (p1, p2) in sort(collect(Tuple{Int,Int}, pr.I))
            v = zeros(ComplexF64, 1 << L); v[1] = 1.0; v[1 + (1 << p1) + (1 << p2)] = 1.0
            itab = mul_tables(itab, v, L)
        end
        ze = Zemul(stab, itab)
        @printf("  K=%-14s I=%-14s : real %10.6f  emul %10.6f  err %8.1e\n",
                pr.K, pr.I, real(zr), real(ze), abs(zr - ze))
    end

    # ---------------- B2：真骨架 × 真IF ----------------
    println("\n── B2：真骨架(εd=0)×真IF(Δ=D_full) ──")
    mpsK0 = GTEMPO.boundarycondition!(GTEMPO.vacuumstate(ComplexF64, lat), lat; band=1, trunc=NoTruncation())
    Sext = mps_to_dense(mpsK0, L)
    Δfull, ν_a = ring_delta(cells)
    mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(Δfull); band=1, trunc=NoTruncation())
    Iext = mps_to_dense(mpsI, L)
    zr = GTEMPO.integrate(lat, mpsK0, mpsI)
    ze = Zemul(Sext, Iext)
    @printf("  real %10.6f  emul %10.6f  err %8.1e\n", real(zr), real(ze), abs(zr - ze))

    # ---------------- B3：εd=0.7 门情形（gate_sign 标定）----------------
    println("\n── B3：εd=0.7 解析 K 链 × 真IF（gate_sign 标定 vs 真管线）──")
    zreal = run_pipeline_real(Δfull, 1.0, EPSD; N=N)
    for gs in (+1.0, -1.0)
        Sk = analytic_kchain(EPSD, nl, lat, L; gate_sign=gs)
        ze = Zemul(Sk, Iext)
        @printf("  gate_sign=%+.0f : 真管线 %10.6f  τ-模拟 %10.6f  err %8.1e\n",
                gs, real(zreal), real(ze), abs(zreal - ze))
    end
    # ---------------- B3'：真 K 链的 h-向量读出（τ2 泛函的 I 侧边际）----------------
    """h(b) := Z_real(K_real(εd), 单项式 b)。g(B)=Z(K, inj(B)) 的 I 侧 Möbius 读全。
    恒等式：Z(K, I) = Σ_b I_b·h(b)（τ2 双线性，不需解 K 系数向量）。"""
    function h_readout(eps_d; N=2, order=2)
        sched = NativeExact.native_schedule(order, N, BETA / N)
        ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
        nlr = length(ring_links)
        latr = GTEMPO.GrassmannLattice(N=nlr, δτ=BETA / nlr, contour=:imag, bands=1)
        mpsK = GTEMPO.vacuumstate(ComplexF64, latr)
        model = ToulouseIM(ϵ_d=eps_d)
        for g in 1:nlr
            fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
            sparse = GTEMPO._propagator_sparsegmps(latr, fm, g, :τ, nothing, nothing)
            mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
        end
        mpsK = GTEMPO.boundarycondition!(mpsK, latr; band=1, trunc=NoTruncation())
        gB = [GTEMPO.integrate(latr, mpsK, injected_mps(latr, pairs_of(B))) for B in evensets]
        h = zeros(ComplexF64, ns)
        for (ib, B) in enumerate(evensets)
            val = 0.0 + 0.0im
            for Bsub in all_subsets(pairs_of(B))
                sb = sort(collect(Iterators.flatten(Bsub)))
                sb = isempty(sb) ? Int[] : collect(sb)
                sgn = (-1.0)^((length(B) - length(sb)) ÷ 2)
                val += sgn * gB[idxof(sb)]
            end
            h[ib] = val
        end
        return h, latr, ring_links, mpsK
    end

    Δfull, ν_a = ring_delta(cells)
    mpsI_D = GTEMPO.hybriddynamics(lat, ImagDelta(Δfull); band=1, trunc=NoTruncation())
    Iext = mps_to_dense(mpsI_D, L)

    """Z(εd, Δ) = Σ_b I(Δ)_b·h_{εd}(b)（精确模拟器）。"""
    function z_em(hvec, Δ)
        IΔ = analytic_if(Δ, lat, L)
        s = 0.0 + 0.0im
        for idx in 1:length(IΔ)
            IΔ[idx] == 0 && continue
            b = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
            ib = findfirst(==(b), evensets)
            ib === nothing && continue
            s += IΔ[idx] * hvec[ib]
        end
        return s
    end
    function dz_em(hvec, Δ, i, j)
        # ∂I/∂Δ_ij = exp(ā_i a_j)·∏_{r≠i} exp(row_r)：用修改 Δ 的积规则
        Nτ = lat.Nτ
        # 行因子表
        rowtab(r, E) = begin
            monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
            for c in 1:Nτ
                iszero(E[r, c]) && continue
                pb = lat[GrassmannContourIndex(r + 1, :τ, 1, true)] - 1
                pe = lat[GrassmannContourIndex(c + 1, :τ, 1, false)] - 1
                sgn = (pb > pe) ? -1.0 : 1.0
                push!(monos, (sort([pb, pe]), E[r, c] * sgn))
            end
            t = zeros(ComplexF64, 1 << L); t[1] = 1.0
            mul_monomials!(zeros(ComplexF64, 1 << L), t, L, Monomials(monos))
        end
        acc = rowtab(1, Δ)
        for r in 2:Nτ
            acc = mul_tables(acc, rowtab(r, Δ), L)
        end
        # 乘 ∂row_i/∂Δ_ij = 单对 exp
        pb = lat[GrassmannContourIndex(i + 1, :τ, 1, true)] - 1
        pe = lat[GrassmannContourIndex(j + 1, :τ, 1, false)] - 1
        sgn = (pb > pe) ? -1.0 : 1.0
        E1 = zeros(ComplexF64, 1 << L); E1[1] = 1.0; E1[1 + (1 << pb) + (1 << pe)] = sgn
        acc = mul_tables(acc, E1, L)
        s = 0.0 + 0.0im
        for idx in 1:length(acc)
            acc[idx] == 0 && continue
            b = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
            ib = findfirst(==(b), evensets)
            ib === nothing && continue
            s += acc[idx] * hvec[ib]
        end
        return s
    end

    println("\n── B3'：h-模拟器验证（Δ=D_full，εd ∈ {0, 0.7}，order 2）──")
    for ed in (0.0, EPSD)
        h2, latr, rl, mpsKr = h_readout(ed; N=N, order=2)
        zreal = GTEMPO.integrate(latr, mpsKr, mpsI_D)
        zemu = z_em(h2, Δfull)
        @printf("  εd=%.2f : real %10.6f  emul %10.6f  err %8.1e\n", ed, real(zreal), real(zemu), abs(zreal - zemu))
        h1, latr1, rl1, mpsKr1 = h_readout(ed; N=N, order=1)
        sched1 = NativeExact.native_schedule(1, N, BETA / N)
        rl1v = vcat(sched1.links[1] + sched1.links[end], sched1.links[2:end-1])
        cells1 = fill(BETA / N, N)
        Δf1, _ = ring_delta(cells1)   # 同 cells（order 1 的 cells 与 order 2 相同：均匀 N）
        h1v = h1
        # order 1 的锚 Z：
        FZ1 = begin
            m = length(VB); nt = 1 + m; dim = 1 << nt
            ops = [annih_op(dim, k) for k in 0:(nt-1)]
            Hq = zeros(ComplexF64, dim, dim); h = hcoupled()
            for i in 1:nt, j in 1:nt
                iszero(h[i, j]) || (Hq .+= h[i, j] .* (ops[i]' * ops[j]))
            end
            Him = ed .* (ops[1]' * ops[1])
            Zc = Matrix{ComplexF64}(I, dim, dim)
            for p in length(cells1):-1:1
                W = exp(-cells1[p] .* Hq) * exp(-rl1v[p] .* Him)
                Zc = W * Zc
            end
            real(tr(Zc))
        end
        mpsI_D1 = GTEMPO.hybriddynamics(latr1, ImagDelta(Δf1); band=1, trunc=NoTruncation())
        zreal1 = GTEMPO.integrate(latr1, mpsKr1, mpsI_D1)
        zemu1 = z_em(h1v, Δf1)
        @printf("  order1 εd=%.2f : real %10.6f  emul %10.6f  err %8.1e   (Fock 锚 %10.6f)\n",
                ed, real(zreal1), real(zemu1), abs(zreal1 - zemu1), FZ1)
    end


    # ---------------- 落盘：τ2 与 h 向量（供轻进程拟合脚本消费）----------------
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    open(joinpath(outdir, "imag_tau2.csv"), "w") do io
        println(io, "ia,ib,re,im")
        for (ia, A) in enumerate(evensets), (ib, B) in enumerate(evensets)
            v = τ2[(A, B)]
            println(io, "$ia,$ib,$(real(v)),$(imag(v))")
        end
    end
    open(joinpath(outdir, "imag_evensets.csv"), "w") do io
        println(io, "id,set")
        for (i, A) in enumerate(evensets)
            println(io, "$i,$(join(A, ' '))")
        end
    end
    println("τ2 已写入 imag_tau2.csv / imag_evensets.csv")
    flush(stdout)
end
main()
