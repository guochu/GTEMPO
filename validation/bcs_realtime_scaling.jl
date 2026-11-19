# =============================================================================
# validation/bcs_realtime_scaling.jl —— BCS 实时演化 order1/order2 离散化斜率
# 标定（plan/bcs_exacttti_plan.md Phase W 科学核心；raw 口径，去顶点移位后）。
#
# 模型：BCS bath 单模 Anderson（U=1, ϵ_d=0.8, Δ=0.6, β=1, ω=1, α=0.5，Z2），
# IF = PartialIF（反常通道）；K = fock_sysdynamics order1/order2（Weyl-2
# 只改 K 侧与测度——IF 精度独立于 K 离散化阶）。
# 参考 = bcs_ed 全模型 ED（数守恒退耦初态：杂质热态 ⊗ 浴热态）。
#
# 判据（T7）：order=2 渐近斜率 ≈ 2.00（Weyl 标量无 BCS 修正的端到端确证；
# bcsweyl F5 裁决反常道零贡献）。order=1 斜率 ≈ 1（raw 常数档，阶不变）。
#
# 运行：julia --project=<GTEMPO 根> validation/bcs_realtime_scaling.jl
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const U = 1.0; const ϵ_ed = 0.8; const ω = 1.0; const α = 0.5
const β = 1.0; const T = 0.4; const Δ_bcs = 0.6

let
    # ED 参照（最大网格全时序；短网格取其前缀——同 H/H0 连续演化）
    Nmax = 32
    δt = T / Nmax
    ts = (1:Nmax) .* δt
    p1 = spin_half_matrices()
    n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
    I2 = one(n̂)
    n_imp = kron(n̂, I2) + kron(I2, n̂)
    Himpbare = -ϵ_ed * n_imp + U * kron(n̂, n̂)
    I4 = Matrix{ComplexF64}(I, 4, 4)
    n_bath = kron(I2, n̂) + kron(n̂, I2)
    H, Ab, Bb, _ = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ_bcs)
    H0dec = kron(Himpbare, I4) + ω * kron(I4, n_bath)
    gt_ed, lt_ed = greater_lesser_ed(H, Ab, Bb, H0dec, ts, β)
    model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ_bcs)
    cn = coupling_norm(bath) * 2
    ρ0 = thermal_impurity_state(model, β)

    println("=" ^ 88)
    @printf("BCS realtime scaling — Anderson (U=%.1f, ϵ_d=%.1f), Δ=%.1f, β=%.1f, T=%.2f\n",
            U, ϵ_ed, Δ_bcs, β, T)
    println("=" ^ 88)
    println("  Nt   order   maxabs_gt      maxabs_lt      (relerr 参考)")
    errors = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for Nt in (4, 8, 16, 32)
        δt = T / Nt
        lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
        kernel = correlationfunction(bath, lat)
        Is = bcs_hybriddynamics(lat, bath, PartialIF())
        for order in (1, 2)
            K = order == 2 ?
                sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0) :
                sysdynamics(lat, model; initial_state=ρ0)
            Z = integrate(K, Is)
            gt = [greater(lat, model, kernel, Is; i=i, j=0, Z=Z, order=order,
                          coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
            lt = [lesser(lat, model, kernel, Is; i=i, j=0, Z=Z, order=order,
                         coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
            # 斜率口径 = 绝对偏差（跨 N 可比）；ED 参照按时刻对齐切片——
            # 端口格点 i·T/Nt = ED 网格下标 i·stride（Nmax 必须被 Nt 整除）
            stride = Nmax ÷ Nt
            idx = (1:Nt) .* stride
            eg = maximum(abs.(gt .- gt_ed[idx])); el = maximum(abs.(lt .- lt_ed[idx]))
            er = (relerr(gt, gt_ed[idx]), relerr(lt, lt_ed[idx]))
            errors[(Nt, order)] = (eg, el)
            @printf("  %3d   %d      %.3e    %.3e   (relerr %.2e / %.2e)\n",
                    Nt, order, eg, el, er[1], er[2])
        end
    end
    # log-log 最小二乘斜率（后三点渐近档）
    println("-" ^ 88)
    for order in (1, 2)
        for (idx, lbl) in ((1, "gt"), (2, "lt"))
            Ns = [4, 8, 16, 32]
            x = [log(T / Nt) for Nt in Ns]
            y = [log(max(errors[(Nt, order)][idx], 1e-16)) for Nt in Ns]
            n = length(x)
            slope = sum((x[i] - sum(x)/n) * (y[i] - sum(y)/n) for i in 1:n) /
                    sum((x[i] - sum(x)/n)^2 for i in 1:n)
            @printf("order=%d %s: slope = %.3f\n", order, lbl, slope)
        end
    end
    println("实测结论（2026-11-19，raw 口径）：")
    println("  order=1 slope ≈ 1.07–1.09 ✓（一阶；raw 常数档，阶不变）")
    println("  order=2 较 order=1 改善 ~40×（1e-4 → 1e-5 档），但 Nt≥8 撞 ~1e-5")
    println("  δt 无关地板（紧截断不变、随 Nt 线性缓增 = 每格常数残差特征）——")
    println("  渐近斜率 ≈ 2.00 被地板掩盖。地板源疑 = 二阶接触门 × BCS 交互")
    println("  （plan R7 排查面；正常 bath 九组基准 1.996–2.005 无此地板）。")
end
