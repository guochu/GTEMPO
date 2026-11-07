# =============================================================================
# imag_calib_solve.jl —— 消费 imag_trace_oracle.jl 落盘的 τ2/h 数据做拟合
#
# 轻进程：无 GTEMPO 收缩，只有稠密多项式代数（analytic_if / mul_monomials!）。
# 阶段：
#   C   拟合 (Δ, ν) vs Fock 锚（LM，解析 Jacobian）
#   C2  K 链 + 跨位分辨率因子（s = ±1 旋钮）× 拟合 Δ（h 数据须含 C2 变体行）
#   D   结论汇总
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_calib_solve.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO                    # formb_dense 依赖
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact
using QuAPI: AbstractCorrelationFunction

# ------------------------------------------------------------------ 模型 --
const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = begin
    m = length(VB)
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = conj(VB[k]); h[k+1, 1] = VB[k]; h[k+1, k+1] = HB[k, k]
    end
    h
end
annih_op(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end
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
function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim); h = hcoupled()
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
    return real(tr(Z)), real(tr(Zn) / tr(Z))
end

# ------------------------------------------------------- ImagDelta 适配器 --
mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
end
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex) = x.Δ[i.j, j.j]

"""真管线（可含 K 侧跨位分辨率因子 s·pairs）。"""
function run_pipeline_real(Δ, ν, eps_d; N=2, order=2, res_s=0.0)
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
    if res_s != 0
        for p in 1:N
            s1 = p + 1
            s2 = p == N ? 2 : p + 2
            pb = lat[GrassmannContourIndex(s1, :τ, 1, true)]
            pe = lat[GrassmannContourIndex(s2, :τ, 1, false)]
            mpsK = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(res_s))), mpsK)
        end
    end
    Z = GTEMPO.integrate(lat, mpsK, mpsI) * ν
    b1 = GrassmannContourIndex(1, conj=true, branch=:τ, band=1)
    g_end = real(GTEMPO.gf(lat, (GrassmannContourIndex(lat.Nτ, conj=false, branch=:τ, band=1), b1),
                           mpsK, mpsI, Z=Z))
    g_1 = real(GTEMPO.gf(lat, (GrassmannContourIndex(1, conj=false, branch=:τ, band=1), b1),
                         mpsK, mpsI, Z=Z))
    return Z, (g=g_end, g1=g_1, om1=1 - g_1, ome=1 - g_end)
end

# ------------------------------------------------------------------ 数据 --
const NL = 2
const N = 2
const CELLS = fill(BETA / N, N)
const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
const L = length(LAT)
const EVENS = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const NS = length(EVENS)
const IDXOF = Dict(A => i for (i, A) in enumerate(EVENS))
const GRID = collect(0.0:0.125:1.875)
const HOLDOUT = [0.525, 1.225]     # 0.525 在网格内；留出 = 用 0.5625/1.1875 则需另跑
const HOLDOUT2 = [0.0625, 1.9375]  # （仅占位——实际留出检验见 D 阶段说明）

function load_data()
    dir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    τ2 = Dict{Tuple{Int,Int}, ComplexF64}()
    for (i, ln) in enumerate(eachline(joinpath(dir, "imag_tau2.csv")))
        i == 1 && continue
        f = split(strip(ln), ',')
        τ2[(parse(Int, f[1]), parse(Int, f[2]))] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
    end
    hdata = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()   # (order, epsd) → h
    for (i, ln) in enumerate(eachline(joinpath(dir, "imag_h.csv")))
        i == 1 && continue
        f = split(strip(ln), ',')
        key = (parse(Int, f[1]), parse(Float64, f[2]))
        ib = parse(Int, f[3])
        get!(hdata, key, zeros(ComplexF64, NS))[ib] = ComplexF64(parse(Float64, f[4]), parse(Float64, f[5]))
    end
    return τ2, hdata
end

# ------------------------------------------- 稠密 IF 表与解析 Jacobian --
const BARPOS(i) = LAT[GrassmannContourIndex(i + 1, :τ, 1, true)] - 1
const ETAPOS(j) = LAT[GrassmannContourIndex(j + 1, :τ, 1, false)] - 1
const PAIRSGN(i, j) = BARPOS(i) > ETAPOS(j) ? -1.0 : 1.0
const ROWPOS = [(BARPOS(i), ETAPOS(j), PAIRSGN(i, j)) for i in 1:N, j in 1:N]

"""I(Δ) 表：∏ᵢ exp(ā_i Σ_j Δ_ij a_j)（canonical）。"""
function if_table(Δ)
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for i in 1:N
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        for j in 1:N
            iszero(Δ[i, j]) && continue
            pb, pe, sg = ROWPOS[i, j]
            push!(monos, (sort([pb, pe]), Δ[i, j] * sg))
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end

"""∂I/∂Δ_ij = 单对 exp × ∏_{r≠i} exp(row_r)。"""
function dif_table(Δ, i, j)
    acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
    for r in 1:N
        monos = Tuple{Vector{Int}, ComplexF64}[(Int[], 1.0 + 0.0im)]
        if r == i
            pb, pe, sg = ROWPOS[i, j]
            push!(monos, (sort([pb, pe]), sg))
        else
            for c in 1:N
                iszero(Δ[r, c]) && continue
                pb, pe, sg = ROWPOS[r, c]
                push!(monos, (sort([pb, pe]), Δ[r, c] * sg))
            end
        end
        acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L, Monomials(monos))
    end
    return acc
end

"""canonical 位集多项式乘（表 → 表）。"""
function mul_tables(a::Vector{ComplexF64}, b::Vector{ComplexF64}, L::Int)
    monos = Tuple{Vector{Int}, ComplexF64}[]
    for idx in 1:length(b)
        b[idx] == 0 && continue
        bits = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
        push!(monos, (bits, b[idx]))
    end
    return mul_monomials!(zeros(ComplexF64, 1 << L), a, L, Monomials(monos))
end

"""Z = Σ_b I_b·h_b（h 键 = 偶集合索引）。"""
function z_em(hv, Itab)
    s = 0.0 + 0.0im
    for idx in 1:length(Itab)
        Itab[idx] == 0 && continue
        b = [p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
        ib = get(IDXOF, b, 0)
        ib == 0 && continue
        s += Itab[idx] * hv[ib]
    end
    return s
end

# ================================================================== 主程序 --
function main()
    τ2, hdata = load_data()
    nzτ = count(v -> abs(v) > 1e-12, values(τ2))
    @printf("τ2：%d 条目（非零 %d）\n", length(τ2), nzτ)

    Δfull, ν_a = ring_delta(CELLS)
    println("锚（Fock 环）：")
    zanch = Dict{Float64, Float64}()
    for ed in vcat(GRID, HOLDOUT)
        cells2 = CELLS; durs2 = CELLS   # order 2 均匀：cells=links=δ
        z, n = fock_anchor(cells2, durs2, ed)
        zanch[ed] = z
        @printf("  εd=%.3f : Z = %.12f  n = %.12f\n", ed, z, n)
    end
    @printf("D_full = %s\nν_a = %.10f\n", round.(Δfull, digits=6), real(ν_a))

    # ---------------- C：拟合 (Δ, ν) ----------------
    println("\n── C：拟合 x = [Δ11,Δ12,Δ21,Δ22,ν] ──")
    function resid(x)
        Δ = ComplexF64.(reshape(x[1:4], 2, 2))
        return [x[5] * real(z_em(hdata[(2, ed)], if_table(Δ))) - zanch[ed] for ed in GRID]
    end
    function jacob(x)
        Δ = ComplexF64.(reshape(x[1:4], 2, 2))
        J = zeros(length(GRID), 5)
        for (r, ed) in enumerate(GRID)
            hv = hdata[(2, ed)]
            for (p, (i, j)) in enumerate([(1, 1), (1, 2), (2, 1), (2, 2)])
                J[r, p] = x[5] * real(z_em(hv, dif_table(Δ, i, j)))
            end
            J[r, 5] = real(z_em(hv, if_table(Δ)))
        end
        J
    end
    x = [real.(Δfull)...; 1.0]
    λ = 1e-6
    for it in 1:200
        r = resid(x); J = jacob(x)
        dx = (J' * J + λ * diagm(0 => max.(diag(J' * J), 1e-12))) \ (-J' * r)
        x2 = x + dx
        if norm(resid(x2)) < norm(r)
            x = x2; λ = max(λ / 3, 1e-14)
        else
            λ *= 5
        end
        it % 10 == 0 && @printf("  iter %3d : |r| = %9.3e  x = %s\n", it, norm(r), [round(v, digits=5) for v in x])
        norm(r) < 1e-12 && break
    end
    Δc = ComplexF64.(reshape(x[1:4], 2, 2))
    @printf("C 结果：|r| = %.3e  ν = %.8f\nΔ = ", norm(resid(x)), x[5])
    display(round.(Δc, digits=8))
    rh = [x[5] * real(z_em(hdata[(2, ed)], if_table(Δc))) - zanch[ed] for ed in HOLDOUT]
    @printf("留出残差 = %s\n", [round(abs(v), sigdigits=3) for v in rh])
    println("（对照 ν_a = ", round(real(ν_a), digits=6), "）")
    flush(stdout)

    # ---------------- C2：I 侧 = row-exp(Δ) × 跨位分辨率因子 ----------------
    println("\n── C2：I(Δ) × 分辨率因子（kind ∈ {bar→eta, eta→bar} × s ∈ {±1}）──")
    """跨位分辨率对表：p=1..N，kind=1: pair(bar_{slice p+1}, eta_{slice p+2})；
    kind=2: pair(eta_{slice p+1}, bar_{slice p+2})；环闭合（p=N → 对方切片 2）。"""
    function res_table(s::Real, kind::Int)
        acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
        for p in 1:N
            s1 = p + 1
            s2 = p == N ? 2 : p + 2
            if kind == 1
                pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)] - 1
                pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)] - 1
            else
                pb = LAT[GrassmannContourIndex(s2, :τ, 1, true)] - 1
                pe = LAT[GrassmannContourIndex(s1, :τ, 1, false)] - 1
            end
            sg = (pb > pe) ? -1.0 : 1.0
            acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L,
                                 Monomials([(Int[], 1.0 + 0.0im),
                                            (sort([pb, pe]), ComplexF64(s * sg))]))
        end
        return acc
    end
    function lm_fit(fit_resid, fit_jacob, x0; iters=200)
        y = copy(x0); λ = 1e-6
        for it in 1:iters
            r = fit_resid(y); J = fit_jacob(y)
            dx = (J' * J + λ * diagm(0 => max.(diag(J' * J), 1e-12))) \ (-J' * r)
            y2 = y + dx
            if norm(fit_resid(y2)) < norm(r)
                y = y2; λ = max(λ / 3, 1e-14)
            else
                λ *= 5
            end
            norm(r) < 1e-13 && break
        end
        return y, norm(fit_resid(y))
    end
    best = nothing
    for kind in (1, 2), s in (1.0, -1.0)
        rt = res_table(s, kind)
        function resid2(x)
            Δ = ComplexF64.(reshape(x[1:4], 2, 2))
            itab = mul_tables(if_table(Δ), rt, L)
            return [x[5] * real(z_em(hdata[(2, ed)], itab)) - zanch[ed] for ed in GRID]
        end
        function jacob2(x)
            Δ = ComplexF64.(reshape(x[1:4], 2, 2))
            J = zeros(length(GRID), 5)
            for (r, ed) in enumerate(GRID)
                hv = hdata[(2, ed)]
                for (p, (i, j)) in enumerate([(1, 1), (1, 2), (2, 1), (2, 2)])
                    J[r, p] = x[5] * real(z_em(hv, mul_tables(dif_table(Δ, i, j), rt, L)))
                end
                itab = mul_tables(if_table(Δ), rt, L)
                J[r, 5] = real(z_em(hv, itab))
            end
            J
        end
        x0 = [real.(Δfull)...; 1.0]
        xf, rf = lm_fit(resid2, jacob2, x0)
        @printf("  kind=%d s=%+.0f : |r| = %9.3e  Δ = %s  ν = %.6f\n", kind, s, rf,
                [round(v, digits=4) for v in xf], xf[5])
        if best === nothing || rf < best[1]
            best = (rf, kind, s, copy(xf))
        end
        # 多起点（D_full×随机扰动）防坏盆地
        for trial in 1:8
            y0 = [real.(Δfull)...] .+ 0.8 .* (rand(4) .- 0.5)
            xf2, rf2 = lm_fit(resid2, jacob2, [y0..., 1.0])
            if rf2 < best[1]
                best = (rf2, kind, s, copy(xf2))
            end
            rf2 < 1e-10 && break
        end
    end
    println("C2 最优：|r| = ", best[1], "  kind=", best[2], "  s=", best[3])
    println("  x = ", [round(v, digits=6) for v in best[4]])

    # ---------------- C3：逐对分辨率系数 c_p 自由化（LM 多起点）----------------
    println("\n── C3：I(Δ) × ∏_p exp(c_p·pair_p)（c_p 自由复数，kind=1 位）──")
    # 参数 x = [Δ11re,Δ11im,Δ12re,Δ12im,Δ21re,Δ21im,Δ22re,Δ22im, c1re,c1im,c2re,c2im, ν]
    function res_pair_table(cs::Vector{ComplexF64})
        acc = zeros(ComplexF64, 1 << L); acc[1] = 1.0
        for p in 1:N
            s1 = p + 1
            s2 = p == N ? 2 : p + 2
            pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)] - 1
            pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)] - 1
            sg = (pb > pe) ? -1.0 : 1.0
            acc = mul_monomials!(zeros(ComplexF64, 1 << L), acc, L,
                                 Monomials([(Int[], 1.0 + 0.0im),
                                            (sort([pb, pe]), cs[p] * sg)]))
        end
        return acc
    end
    function unpack(x)
        Δ = [ComplexF64(x[1], x[2]) ComplexF64(x[3], x[4]);
             ComplexF64(x[5], x[6]) ComplexF64(x[7], x[8])]
        cs = [ComplexF64(x[9], x[10]), ComplexF64(x[11], x[12])]
        return Δ, cs, x[13]
    end
    function itab_of(x)
        Δ, cs, ν = unpack(x)
        return mul_tables(if_table(Δ), res_pair_table(cs), L)
    end
    function resid3(x)
        itab = itab_of(x)
        return [x[13] * real(z_em(hdata[(2, ed)], itab)) - zanch[ed] for ed in GRID]
    end
    function jacob3(x)
        Δ, cs, ν = unpack(x)
        rt = res_pair_table(cs)
        itab = mul_tables(if_table(Δ), rt, L)
        J = zeros(length(GRID), 13)
        dpos = [(1, 1), (1, 2), (2, 1), (2, 2)]
        for (r, ed) in enumerate(GRID)
            hv = hdata[(2, ed)]
            for (p, (i, j)) in enumerate(dpos)
                dz = z_em(hv, mul_tables(dif_table(Δ, i, j), rt, L))
                J[r, 2p-1] = ν * real(dz)
                J[r, 2p] = ν * real(-im * dz)
            end
            for p in 1:N
                s1 = p + 1
                s2 = p == N ? 2 : p + 2
                pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)] - 1
                pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)] - 1
                sg = (pb > pe) ? -1.0 : 1.0
                E1 = zeros(ComplexF64, 1 << L); E1[1] = 1.0; E1[1 + (1 << pb) + (1 << pe)] = sg
                rest = res_pair_table([q == p ? ComplexF64(0) : cs[q] for q in 1:N])
                dtab = mul_tables(mul_tables(if_table(Δ), rest, L), E1, L)
                dz = z_em(hv, dtab)
                J[r, 8 + 2p - 1] = ν * real(dz)
                J[r, 8 + 2p] = ν * real(-im * dz)
            end
            J[r, 13] = real(z_em(hv, itab))
        end
        return J
    end
    found = 0
    global best3 = (Inf, nothing)
    for trial in 1:30
        x0 = vcat(vec(real.(Δfull)), 0.0, 0.0, 0.0, 0.0,
                  1.0 - rand(), 0.2 * (rand() - 0.5), 1.0 - rand(), 0.2 * (rand() - 0.5), 1.0)
        y, rf = lm_fit(resid3, jacob3, x0; iters=400)
        if rf < best3[1]
            global best3 = (rf, copy(y))
        end
        if rf < 1e-10
            global found += 1
            Δ, cs, ν = unpack(y)
            @printf("  trial %d: |r| = %8.1e  Δ = %s  c = %s  ν = %.6f\n", trial, rf,
                    [round(real(v), digits=4) for v in Δ], [round(real(v), digits=4) for v in cs], ν)
            found >= 4 && break
        end
    end
    println("C3 命中（|r|<1e-10）次数 = ", found, "；最优 |r| = ", best3[1])
    if best3[1] < 1e-6
        Δ, cs, ν = unpack(best3[2])
        println("最优解：")
        display(round.(Δ, digits=6))
        println("c = ", [round(v, digits=6) for v in cs], "  ν = ", ν)
        println("对照 D_full = ", round.(Δfull, digits=6))
    end
end
main()
