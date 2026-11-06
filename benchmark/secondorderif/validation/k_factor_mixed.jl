# ============================================================================
# k_factor_mixed.jl — 混合 K 的 U=0 因子化逐支撑楔积判决（第十轮核心工具）
#
# 用法：julia --project=<repo> benchmark/secondorderif/validation/k_factor_mixed.jl <repo>
#
# Kshould = K₁（单带 mixed，已验证 clean）两条带嵌入的 graded 楔积：
#   φ(m₁∪m₂) = φ₁(m₁)·φ₂(m₂)·(−1)^{inv(m)}，
#   inv(m) = #{(x∈m₁, y∈m₂) : pos_x > pos_y}（ChronoMixed 块内带 1 先带 2，
#   同格对永不反转 ⟹ 只统计跨格反转）。
# 与 Q = coeffs(K₂(U=0)) 逐单项比对，失配按 (band1 变量数, band2 变量数) 分类。
#
# 判决史（2026-11-06）：修复前 96/400 单项反号且全部 deg₁≥2∧deg₂≥2——
# 总宇称补偿 (−1)^{c_total·p_total} 相对每带 XOR ⊕_b c_b·p_b(l_b) 恰差交叉项
# c₁p₂⊕c₂p₁（只在两带都选中时非零）。system.jl 改每带 XOR 后 **0/400**。
# ⚠️ 内存：L=24 提取（2^24 ≈ 270MB）安全。
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using GTEMPO
const TK = GTEMPO.TK
using Printf
using LinearAlgebra: mul!

μ, ω, α, β = 0.7, 1.0, 0.5, 1.0

function coeffs(mps)
    L = length(mps)
    Ts = map(1:L) do i
        t = mps[i]
        d0i = TK.blockdim(TK.space(t, 1), TK.Z2Irrep(0)); d1i = TK.blockdim(TK.space(t, 1), TK.Z2Irrep(1))
        d0o = TK.blockdim(TK.domain(t)[1], TK.Z2Irrep(0)); d1o = TK.blockdim(TK.domain(t)[1], TK.Z2Irrep(1))
        T0 = zeros(ComplexF64, d0i + d1i, d0o + d1o); T1 = zeros(ComplexF64, d0i + d1i, d0o + d1o)
        for (f1, f2) in TK.fusiontrees(t)
            ql, p = Int(f1.uncoupled[1].n), Int(f1.uncoupled[2].n)
            qr = Int(f2.uncoupled[1].n)
            B = t[f1, f2]
            offq = ql == 0 ? 0 : d0i; offr = qr == 0 ? 0 : d0o
            for a in 1:size(B, 1), bb in 1:size(B, 3)
                v = B[a, 1, bb]
                iszero(v) && continue
                T = p == 0 ? T0 : T1
                T[offq + a, offr + bb] += v
            end
        end
        (T0, T1)
    end
    env = reshape(ComplexF64[1.0], 1, 1)
    for i in 1:L
        T0, T1 = Ts[i]
        n = size(env, 1)
        ne = zeros(ComplexF64, 2n, size(T1, 2))
        mul!(view(ne, 1:n, :), env, T0); mul!(view(ne, n+1:2n, :), env, T1)
        env = ne
    end
    co = GTEMPO.scaling(mps) .* env[:, 1]
    d = Dict{Int,ComplexF64}()
    for m in 0:(length(co) - 1)
        abs(co[m+1]) > 1e-10 && (d[m] = co[m+1])
    end
    return d
end

bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)

# ---- 单带 mixed K（U=0 ⟹ ToulouseIM；含已标定的补偿）----
c1 = WeylMixedTime(Nτ = 1, dτ = β, Nt = 1, dt = 0.2, bands = 1)
kern1 = hybridization(bath, c1)
K1 = system_factor(c1, ToulouseIM(ϵ_d = μ); order = 2, coupling_norm = α)
P1 = coeffs(K1)

# ---- 双带 U=0 mixed K ----
c2 = WeylMixedTime(Nτ = 1, dτ = β, Nt = 1, dt = 0.2, bands = 2)
kern2 = hybridization(bath, c2)
K2 = system_factor(c2, AndersonIM(ϵ_d = μ, U = 0.0); order = 2, coupling_norm = α)
Q = coeffs(K2)

lattice1 = weyl_lattice(c1); lattice2 = weyl_lattice(c2)

# 单带 MPS 站点（1-based）→ (cell, conj)；端点不映射
inv1 = Dict{Int,Tuple{Int,Bool}}()
for v in 2:lattice1.kτ
    inv1[lattice1[GrassmannContourIndex(v, :τ, 1, false)]] = (v - 1, false)
    inv1[lattice1[GrassmannContourIndex(v, :τ, 1, true)]] = (v - 1, true)
end
for v in 1:lattice1.Nt, br in (:+, :-), cc in (false, true)
    s = lattice1[GrassmannContourIndex(v, br, 1, cc)]
    (br !== :τ && v == lattice1.kt) && continue
    cell = br === :+ ? lattice1.Nτ + v : lattice1.Nτ + lattice1.Nt + (lattice1.Nt - v + 1)
    inv1[s] = (cell, cc)
end

# 单带 (cell, conj) → 双带各带站点位（1-based）
function posmap(band)
    pm = Dict{Tuple{Int,Bool},Int}()
    for v in 2:lattice2.kτ, cc in (false, true)
        pm[(v - 1, cc)] = lattice2[GrassmannContourIndex(v, :τ, band, cc)]
    end
    for v in 1:lattice2.Nt, br in (:+, :-), cc in (false, true)
        (br !== :τ && v == lattice2.kt) && continue
        cell = br === :+ ? lattice2.Nτ + v : lattice2.Nτ + lattice2.Nt + (lattice2.Nt - v + 1)
        pm[(cell, cc)] = lattice2[GrassmannContourIndex(v, br, band, cc)]
    end
    return pm
end
pm1 = posmap(1); pm2 = posmap(2)

# 双带位串 → [(cell, conj)…] per band（端点位（cell 0 / cell > ncells）跳过 ⟹ 失败标记）
function bits2desc(m)
    d1 = Tuple{Int,Bool}[]; d2 = Tuple{Int,Bool}[]
    mm = m
    while mm != 0
        b = mm & (-mm); mm -= b
        p = trailing_zeros(b) + 1
        if haskey(sm2a_band1, p)
            push!(d1, sm2a_band1[p])
        elseif haskey(sm2a_band2, p)
            push!(d2, sm2a_band2[p])
        else
            return nothing, nothing
        end
    end
    return d1, d2
end
sm2a_band1 = Dict{Int,Tuple{Int,Bool}}()
sm2a_band2 = Dict{Int,Tuple{Int,Bool}}()
for (k, s) in pm1; sm2a_band1[s] = k; end
for (k, s) in pm2; sm2a_band2[s] = k; end

P1map = Dict{Vector{Tuple{Int,Bool}},ComplexF64}()
for (m1, v1) in P1
    desc = Tuple{Int,Bool}[]
    mm = m1; ok = true
    while mm != 0
        b = mm & (-mm); mm -= b
        s0 = trailing_zeros(b) + 1
        haskey(inv1, s0) || (ok = false; break)
        push!(desc, inv1[s0])
    end
    ok && (P1map[sort(desc)] = v1)
end

"K₁∧K₁ 系数：φ₁·φ₂·(−1)^{跨格 band1→band2 反转数}"
function wedge_coeff(d1, d2)
    φ1 = get(P1map, sort(d1), 0.0im)
    φ2 = get(P1map, sort(d2), 0.0im)
    (iszero(φ1) || iszero(φ2)) && return 0.0im
    inv = 0
    for (c1, _) in d1, (c2, _) in d2
        # band2 变量在块内后位 ⟹ 跨格时 pos_y > pos_x 除非 c2 < c1
        (c2 < c1) && (inv += 1)
    end
    return (isodd(inv) ? -1 : 1) * φ1 * φ2
end

function compare()
    nmm = 0; ncmp = 0
    patterns = Dict{String,Int}()
    worst = Float64[]
    for (m, vq) in Q
        d1, d2 = bits2desc(m)
        d1 === nothing && continue
        vs = wedge_coeff(d1, d2)
        (abs(vs) > 1e-12 || abs(vq) > 1e-12) || continue
        ncmp += 1
        if abs(vs - vq) > 1e-9 * max(1.0, abs(vs), abs(vq))
            nmm += 1
            key = "deg1=$(length(d1)),deg2=$(length(d2))"
            patterns[key] = get(patterns, key, 0) + 1
            push!(worst, abs(vs - vq) / max(abs(vs), abs(vq)))
            if nmm <= 12
                @printf("失配 #%-3d %s: K₂=%+.8f%+.8fi  K₁∧K₁=%+.8f%+.8fi  cells1=%s cells2=%s\n",
                        nmm, key, real(vq), imag(vq), real(vs), imag(vs),
                        string(sort(unique([c for (c, _) in d1]))),
                        string(sort(unique([c for (c, _) in d2]))))
            end
        end
    end
    @printf("\n比对单项 %d；失配 = %d；最大相对失配 = %.3e\n失配形态分布：\n", ncmp, nmm,
            nmm > 0 ? maximum(worst) : 0.0)
    for (k, cnt) in sort(collect(patterns); by = x -> -x[2])
        println("  ", k, " × ", cnt)
    end
end
compare()
