# =============================================================================
# BCS 连续 bath BdG 锚标定（test/bcs/continuousbath 容差定标）
# 锚 = 谱中点离散化（M 胞，v² = ∫cell J）+ 全二次 BdG（2 + 2M 模）：
#   初始 = 杂质热态 ⊗ 配对浴热态（F0 = 2(I+e^{β𝓗₀})⁻¹ 自动含浴反常密度）
#   gt(t) = −i·[e^{−i𝓗t}·F0]₁₁（左上 n×n 块），lt(t) = +i·[F0_B·e^{−i𝓗t}]₁₁
# 自检：M=1（单胞中点）应 ≡ bcs_ed 单模演化；M→M/2 锚漂移须 ≪ 测试容差。
# =============================================================================
using GTEMPO
using LinearAlgebra: I, dot
using Printf
using QuadGK: quadgk

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

"中点离散化 BdG 锚：返回 (gt(t), lt(t)) 列表"
function bdg_anchor(spec, Δ; n_imp=2, ϵ_d, β, ts, M, δw=nothing)
    lb, ub = spec.lb, spec.ub
    edges = range(lb, ub; length=M+1)
    ws = [(edges[k]+edges[k+1])/2 for k in 1:M]
    v2 = [quadgk(spec.f, edges[k], edges[k+1]; rtol=1e-10)[1] for k in 1:M]
    nb = length(ws)
    n = n_imp + 2nb                       # 杂质 n_imp + 每浴模 Nambu 对 2 模
    h = zeros(ComplexF64, n, n)
    for i in 1:n_imp; h[i, i] = ϵ_d; end
    Π = zeros(ComplexF64, n, n)
    for k in 1:nb
        h[n_imp+2k-1, n_imp+2k-1] = ws[k]
        h[n_imp+2k, n_imp+2k] = ws[k]
        Π[n_imp+2k-1, n_imp+2k] = -Δ
        Π[n_imp+2k, n_imp+2k-1] = +Δ
        v = sqrt(v2[k])
        h[1, n_imp+2k-1] = h[n_imp+2k-1, 1] = v   # d↑ ↔ b↑_k
        h[2, n_imp+2k] = h[n_imp+2k, 2] = v       # d↓ ↔ b↓_k
    end
    # 退耦 BdG（热态密度用）与耦合 BdG（演化用）
    h0 = copy(h)
    for i in 1:n_imp, k in 1:(n - n_imp)÷2
        h0[i, n_imp+2k-1] = h0[n_imp+2k-1, i] = 0
        h0[i, n_imp+2k] = h0[n_imp+2k, i] = 0
    end
    H0dec = [h0 Π; Π -conj(h0)]
    H0cpl = [h Π; Π -conj(h)]
    # BdG 热态密度（⟨Ξ†Ξ⟩ = (I+e^{β𝓗₀ᵈᵉᶜ})⁻¹ 约定；配对浴反常密度自动含于跨块）
    G = inv(I + exp(β * H0dec))
    Gt = transpose(G)
    out = Vector{Tuple{ComplexF64,ComplexF64}}()
    for t in ts
        Mx = exp(-im * H0cpl * t)         # 耦合 BdG 演化
        # G⁺(t)₁₁ = [M − M·Gᵀ]₁₁；G⁻(t)₁₁ = Σ_k M[1,k]G[1,k] − Σ_k M[1,n+k]·conj(G[1,n+k])
        gp11 = (Mx - Mx * Gt)[1, 1]
        gm11 = dot(Mx[1, 1:n], Vector(G[1, 1:n])) -
               dot(Mx[1, n+1:2n], conj.(Vector(G[1, n+1:2n])))
        push!(out, (-im * gp11, im * gm11))
    end
    return out
end

const J(ϵ) = 0.25 * cos(π * (ϵ - 0.2) / 4)^2      # toulouse 基准同款紧凑谱
spec = spectrum(J, lb=-1.8, ub=2.2)
const BETA = 1.0; const EPSD = 0.7; const DELTA = 0.6; const T = 0.4
ts = (1:16) .* 0.025

# 自检 2：锚 vs bcs_ed——直接构造 ω=1 单模情形（谱 = 窄 Dirac 近似不可行，
# 改用独立途径：M 胞锚在 Δ=0 时退化为 normal 中点锚，与 greater_lesser_ed
# 的单模版（bathspecs=[(ω,α)] 单模）对拍需 δw→0——此处直接对 Δ=0.6 双带
# BdG 做逐块自检：配对浴热态的 ⟨b↓b↑⟩ 解析式 vs F0 反常块）
let Δ = DELTA; β = BETA
    E = sqrt(1.0 + Δ^2)
    κ_exact = -Δ / (2E) * tanh(β * E / 2)      # ⟨b↓ b↑⟩ 解析
    # M=1 锚（ω=1, v²=α=0.5）的 F0 反常块
    n = 4
    h0 = zeros(ComplexF64, n, n)
    h0[1,1] = h0[2,2] = EPSD; h0[3,3] = h0[4,4] = 1.0
    h0[1,3] = h0[3,1] = sqrt(0.5); h0[2,4] = h0[4,2] = sqrt(0.5)
    Π = zeros(ComplexF64, n, n); Π[3,4] = -Δ; Π[4,3] = +Δ
    H0 = [h0 Π; Π -conj(h0)]
    G = inv(I + exp(β * H0))
    κ_G = conj(G[3, 4])                          # ⟨c†↑ c†↓⟩ = conj(G[3,4])（G[ba] = ⟨Ξ_b†Ξ_a⟩）
    println("自检 2：G[3,4] = ", G[3, 4], "  反常密度 = ", κ_G, "  解析 κ = ", κ_exact, "  dev = ", abs(κ_G - κ_exact))
end

# 锚自收敛（M=320 vs 160 vs 80）
println("\n锚自收敛（M=80/160/320，Δ=$DELTA）：")
a80 = bdg_anchor(spec, DELTA; ϵ_d=EPSD, β=BETA, ts=ts, M=80)
a160 = bdg_anchor(spec, DELTA; ϵ_d=EPSD, β=BETA, ts=ts, M=160)
a320 = bdg_anchor(spec, DELTA; ϵ_d=EPSD, β=BETA, ts=ts, M=320)
for (i, t) in enumerate(ts)
    @printf("  t=%.2f  |a80−a160|=%.2e  |a160−a320|=%.2e  gt=%+.6f\n", t,
            abs(a80[i][1] - a160[i][1]), abs(a160[i][1] - a320[i][1]), real(a320[i][1]))
end

# TN（order=2）vs 锚
println("\nTN order=2 vs 锚（M=160）：")
bath = bcsbath(fermionicbath(spec, β=BETA), Δ=DELTA)
model = AndersonIM(U=0.0, ϵ_d=EPSD)
ρ0 = thermal_impurity_state(model, BETA)
cn = coupling_norm(bath) * 2
for Nt in (8, 16), D in (128, 512)
    δt = T / Nt
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF(trunc=GTEMPO.Defaults.alg_trunc_if(D=D)))
    Is = Is isa Tuple ? Is : (Is,)
    K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
    Z = integrate(K, Is...)
    stride = 16 ÷ Nt
    idx = (1:Nt) .* stride
    dev = 0.0; devz = abs(Z - 1)
    for i in 1:Nt
        g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2,
                    coupling_norm=cn, initial_state=ρ0)
        dev = max(dev, abs(g - a320[idx[i]][1]))
    end
    println("  Nt=$Nt D=$D: max|gt−锚| = $(@sprintf("%.3e", dev))   |Z−1| = $(@sprintf("%.2e", devz))  bond=$(maximum(bonddims(Is[1])))")
end

# ===== Δ=0 对照（配对扇区失配定位）=====
println("\n===== Δ=0 对照 =====")
let Δ = 0.0
    a = bdg_anchor(spec, Δ; ϵ_d=EPSD, β=BETA, ts=ts, M=160)
    bath = bcsbath(fermionicbath(spec, β=BETA), Δ=Δ)
    model = AndersonIM(U=0.0, ϵ_d=EPSD)
    ρ0 = thermal_impurity_state(model, BETA)
    cn = coupling_norm(bath) * 2
    Nt = 16; δt = T / Nt
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    Is = Is isa Tuple ? Is : (Is,)
    K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
    Z = integrate(K, Is...)
    for i in (1, 2, 4, 8, 16)
        g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2,
                    coupling_norm=cn, initial_state=ρ0)
        @printf("  i=%2d t=%.3f  gt TN %+.6f%+.6fi  锚 %+.6f%+.6fi  dev %.2e\n", i, i*δt,
                real(g), imag(g), real(a[i][1]), imag(a[i][1]), abs(g - a[i][1]))
    end
end
println("\n===== Δ=0.6 逐点对照 =====")
let Δ = 0.6
    a = bdg_anchor(spec, Δ; ϵ_d=EPSD, β=BETA, ts=ts, M=160)
    bath = bcsbath(fermionicbath(spec, β=BETA), Δ=Δ)
    model = AndersonIM(U=0.0, ϵ_d=EPSD)
    ρ0 = thermal_impurity_state(model, BETA)
    cn = coupling_norm(bath) * 2
    Nt = 16; δt = T / Nt
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    Is = Is isa Tuple ? Is : (Is,)
    K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
    Z = integrate(K, Is...)
    for i in (1, 2, 4, 8, 16)
        g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2,
                    coupling_norm=cn, initial_state=ρ0)
        @printf("  i=%2d t=%.3f  gt TN %+.6f%+.6fi  锚 %+.6f%+.6fi  dev %.2e\n", i, i*δt,
                real(g), imag(g), real(a[i][1]), imag(a[i][1]), abs(g - a[i][1]))
    end
end
