# ============================================================================
# if_band_invariance.jl — 双带 IF₁ vs 单带 IF 的带不变性探针（M≥2 混合工件定位）
#
# 用法：julia --project=<repo> benchmark/secondorderif/validation/if_band_invariance.jl <repo>
# 原理：单带 mixed IF 与双带 mixed IF₁（band-1 位嵌入后）逐掩码比较；
#       U=0 双带物理上 = 两个独立单带 ⟹ 两者应一致。
# ⚠️ 内存纪律：全系数提取器构建 2^L 向量——只在 L ≤ 24 的格点上使用
#    （当前 Nτ=Nt=1 → L=24 ≈ 270MB 安全；曾因 L=36 → 2^36≈1.1TB 两次爆内存）。
# 现状（2026-11-05）：双带 IF₁ 相对单带 IF 有 ~2e-4 结构性小偏移（部分单项式、
#    非全局符号），与观测量级 6e-4 有限 N 残差同源——M≥2 混合装配补偿形式
#    （M=1 标定的 paritybits(m)·inpar[l]）疑不完备，下轮入口。
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using GTEMPO
const TK = GTEMPO.TK
using Printf
using LinearAlgebra: mul!

μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)

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

# 单带 mixed IF（12 站点）
c1 = WeylMixedTime(Nτ = 1, dτ = β, Nt = 1, dt = 0.2, bands = 1)
kern1 = hybridization(bath, c1)
lattice1 = weyl_lattice(c1)
IF1 = SecondOrderWeyl.influence_functional(c1, kern1)[1]
P1 = coeffs(IF1)

# 双带 band-1 IF（24 站点）
c2 = WeylMixedTime(Nτ = 1, dτ = β, Nt = 1, dt = 0.2, bands = 2)
kern2 = hybridization(bath, c2)
lattice2 = weyl_lattice(c2)
IF2 = SecondOrderWeyl.influence_functional(c2, kern2)
P2 = coeffs(IF2[1])

# 嵌入映射：单带 band-1 站点 → 双带 band-1 站点（按 (顶点, 支, conj) 恒等；跳过闲置端点）
inv1 = Dict{Int,Tuple{Int,Symbol,Bool}}()
for v in 2:lattice1.kτ
    inv1[lattice1[GrassmannContourIndex(v, :τ, 1, false)]] = (v, :τ, false)
    inv1[lattice1[GrassmannContourIndex(v, :τ, 1, true)]] = (v, :τ, true)
end
for v in 1:lattice1.Nt, br in (:+, :-)
    inv1[lattice1[GrassmannContourIndex(v, br, 1, false)]] = (v, br, false)
    inv1[lattice1[GrassmannContourIndex(v, br, 1, true)]] = (v, br, true)
end
inv2 = Dict{Tuple{Int,Symbol,Bool},Int}()
for v in 2:lattice2.kτ
    inv2[(v, :τ, false)] = lattice2[GrassmannContourIndex(v, :τ, 1, false)]
    inv2[(v, :τ, true)] = lattice2[GrassmannContourIndex(v, :τ, 1, true)]
end
for v in 1:lattice2.Nt, br in (:+, :-)
    inv2[(v, br, false)] = lattice2[GrassmannContourIndex(v, br, 1, false)]
    inv2[(v, br, true)] = lattice2[GrassmannContourIndex(v, br, 1, true)]
end

embed = Dict{Int,Int}()
for (s1, key) in inv1
    (v, br, cc) = key
    (br === :τ && v == 1) && continue
    (br !== :τ && v == lattice1.kt) && continue   # 端点不映射
    haskey(inv2, (v, br, cc)) || continue
    embed[s1 - 1] = inv2[(v, br, cc)] - 1          # 0-based
end

global function compare_bands()
    bad = 0
    for (m, v) in P1
        m2 = 0; ok = true
        mm = m
        while mm != 0
            b = mm & (-mm); mm -= b
            s0 = trailing_zeros(b)
            if !haskey(embed, s0)
                ok = false; break
            end
            m2 |= 1 << embed[s0]
        end
        ok || continue
        v2 = get(P2, m2, 0.0)
        if abs(v2 - v) > 1e-9 * max(1, abs(v))
            bad += 1
            if bad <= 4
                println("  mask ", m, " → ", m2, ": 单带 ", round(real(v), digits=8),
                        " 双带 ", round(real(v2), digits=8))
            end
        end
    end
    println("单带 IF vs 双带 IF₁（band-1 位嵌入）：失配 = ", bad, " / ", length(P1))
end

compare_bands()
