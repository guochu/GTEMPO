# =============================================================================
# R7 排查第五步修正版：
# ① normal U=0（单带）解析 vs singlemode_ed —— 校准解析机器
# ② BCS U=0 Δ=0.6：解析 BdG vs bcs_ed 的 0⁺ 极限与短时行为
# ③ TN（Δ=0.6, order=2）小 t 的 Im gt —— 跟随解析(−0.5i) 还是 ED(−0.31i)
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const ω = 1.0; const α = 0.5; const β = 1.0; const T = 0.2

"U=0 BdG 解析（n 模：前 n_imp 个杂质 + 后 n_b 个浴模，Π 浴内配对）
返回 gt(t) = −i·[e^{−i𝓗t}·(I−f0)]_{d↑,d↑}、lt = +i·[f0·e^{−i𝓗t}]_{d↑,d↑}"
function bdg_gt(; n_imp, n_b, μ_imp, ωm, αc, Δ, ts, β)
    n = n_imp + n_b
    h = zeros(ComplexF64, n, n)
    for i in 1:n_imp; h[i, i] = μ_imp; end
    for b in 1:n_b; h[n_imp+b, n_imp+b] = ωm; end
    for i in 1:n_imp
        h[i, n_imp+i] = h[n_imp+i, i] = sqrt(αc)   # 第 i 杂化第 i 浴模
    end
    Π = zeros(ComplexF64, n, n)
    if n_b >= 2
        Π[n_imp+1, n_imp+2] = -Δ
        Π[n_imp+2, n_imp+1] = +Δ
    end
    Hbdg = [h Π; Π -conj(h)]
    f0 = zeros(ComplexF64, n, n)
    for i in 1:n_imp; f0[i, i] = 1 / (exp(β * μ_imp) + 1); end   # 无化学势平移（ρ ∝ e^{−βH0}）
    for b in 1:n_b; f0[n_imp+b, n_imp+b] = 1 / (exp(β * ωm) + 1); end      # 浴 μ=0
    Utop = [exp(-im * Matrix(Hbdg) * t)[1:n, 1:n] for t in ts]
    gt = [(-im * (U * (I - f0)))[1, 1] for U in Utop]
    lt = [( im * (f0 * U))[1, 1] for U in Utop]
    return gt, lt
end

# ---- ① normal 单带校准 ----
let ts = (1:8) .* 0.025
    μ = 0.7
    gt_an, lt_an = bdg_gt(n_imp=1, n_b=1, μ_imp=μ, ωm=ω, αc=α, Δ=0.0, ts=ts, β=β)
    H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
    gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
    println("① normal U=0 单带：解析 vs singlemode_ed")
    for (i, t) in enumerate(ts)
        @printf("  t=%.3f  gt 解析 %+.6f%+.6fi  ED %+.6f%+.6fi  dev %.2e\n", t,
                real(gt_an[i]), imag(gt_an[i]), real(gt_ed[i]), imag(gt_ed[i]),
                abs(gt_an[i] - gt_ed[i]))
    end
    println("  max dev = ", maximum(abs.(gt_an .- gt_ed)))
end

# ---- ② BCS U=0：解析 vs bcs_ed 0⁺ 极限 ----
let Δ = 0.6; ϵ_d = -0.8
    ts = (1:8) .* 0.025
    gt_an, lt_an = bdg_gt(n_imp=2, n_b=2, μ_imp=ϵ_d, ωm=ω, αc=α, Δ=Δ, ts=ts, β=β)
    H, Ab, Bb, _ = bcs_ed(0.0, -ϵ_d; ω₀=ω, α=α, Δ=Δ)
    n̂ = spin_half_matrices()["n"]; I2 = one(n̂)
    n_imp = kron(n̂, I2) + kron(I2, n̂)
    I4 = Matrix{ComplexF64}(I, 4, 4); n_bath = kron(I2, n̂) + kron(n̂, I2)
    H0dec = kron(ϵ_d * n_imp, I4) + ω * kron(I4, n_bath)
    gt_ed, lt_ed = greater_lesser_ed(H, Ab, Bb, H0dec, ts, β)
    println("\n② BCS U=0 Δ=", Δ, " 解析 vs bcs_ed")
    for (i, t) in enumerate(ts)
        @printf("  t=%.3f  gt 解析 %+.6f%+.6fi  ED %+.6f%+.6fi  dev %.2e\n", t,
                real(gt_an[i]), imag(gt_an[i]), real(gt_ed[i]), imag(gt_ed[i]),
                abs(gt_an[i] - gt_ed[i]))
    end
    println("  max dev = ", maximum(abs.(gt_an .- gt_ed)))
    # ②′ bcs_ed 在 Δ=0 时 vs 解析（Δ=0 退化两条独立 normal 带）
    gt_an0, _ = bdg_gt(n_imp=2, n_b=2, μ_imp=ϵ_d, ωm=ω, αc=α, Δ=0.0, ts=ts, β=β)
    H0b, Ab0, Bb0, _ = bcs_ed(0.0, -ϵ_d; ω₀=ω, α=α, Δ=0.0)
    gt_ed0, lt_ed0 = greater_lesser_ed(H0b, Ab0, Bb0, H0dec, ts, β)
    println("  Δ=0 交叉验证：解析 vs bcs_ed max dev = ", maximum(abs.(gt_an0 .- gt_ed0)))
end

# ---- ③ TN 小 t 的 Im gt（Δ=0.6, order=2）----
let Δ = 0.6; ϵ_d = -0.8; Nt = 16
    δt = T / Nt
    ts = (1:Nt) .* δt
    gt_an, _ = bdg_gt(n_imp=2, n_b=2, μ_imp=ϵ_d, ωm=ω, αc=α, Δ=Δ, ts=ts, β=β)
    model = AndersonIM(U=0.0, ϵ_d=ϵ_d)
    ρ0 = thermal_impurity_state(model, β)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
    cn = coupling_norm(bath) * 2
    lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
    kernel = correlationfunction(bath, lat)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
    Z = integrate(K, Is)
    gt = [greater(lat, model, kernel, Is; i=i, j=0, Z=Z, order=2,
                  coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
    println("\n③ TN Δ=", Δ, " order=2（Nt=", Nt, "）小 t vs 解析：")
    for i in 1:4
        @printf("  t=%.3f  gt TN %+.6f%+.6fi  解析 %+.6f%+.6fi  dev %.2e\n", ts[i],
                real(gt[i]), imag(gt[i]), real(gt_an[i]), imag(gt_an[i]),
                abs(gt[i] - gt_an[i]))
    end
    println("  max dev = ", maximum(abs.(gt .- gt_an)))
end
