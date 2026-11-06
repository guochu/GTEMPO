# ============================================================================
# twopoint_m2_frontier.jl — M≥2 跨格奇插入两点函数的边界判决工具（第十轮收尾）
#
# 判决（2026-11-06）：
#   mixed M=2：α→0（IF≈恒等）greater(t,0) vs ED 机器精度全对；α=0.5 错 ~1e-1
#     （U=0 亦错）⟹ 病在非常数扇区（K 的 slot 单项 × IF 收缩）；
#   rt M=2：α→0 也错 ~3e-1 ⟹ 常数路径即病（fold 键的每带偏移）。
#   消元：补偿形式（每带 XOR/组合）不敏感、oddlater 无关、标签 F 必要
#     （off 更糟）⟹ 病根 = 每带奇插入宇称偏移 × 联合 Z2 单比特键表示的
#     原理性冲突，修复需表示层改动（Z2×Z2 或装配重设计）——下轮入口。
#   ⚠️ coeffs() 全多项式提取器对 F 偏移 MPS 失效（单插入提取全空）——
#     K 层多项式判决目前不可用；α 扫描端到端判决是唯一可靠视野。
# 防守：raw_correlation 对 M≥2 且 i≠j 显式抛错（measurement.jl）。
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: raw_correlation
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using Printf

const μ, ω, β, T = 0.7, 1.0, 1.0, 0.4
const bath0 = fermionicbath(DiracDelta(ω = ω, α = 1e-12), β = β)

H, d, _ = SingleModeED.model(μ = μ, U = 0.8, bathspecs = [(ω, 0.0), (ω, 0.0)])
ρc = SingleModeED.thermal(H, β)

for (label, mk) in (("mixed", (N, dt) -> WeylMixedTime(Nτ = N, dτ = β / N, Nt = N, dt = dt, bands = 2)),
                    ("rt", (N, dt) -> WeylRealTime(Nt = N, dt = dt, bands = 2)))
    println("== $label M=2 α→0 ==")
    N = 4; dt = T / N
    sol = solve(Problem(mk(N, dt), AndersonIM(ϵ_d = μ, U = 0.8); bath = bath0, order = 2))
    for i in (1, 2, 4)
        nΔ, gt, lt = SingleModeED.observables(H, d, [(i) * dt], ρc)
        gΔ = gt[1] / (-1.0im)
        gmine = raw_correlation(sol, i, 0; kind = :greater)
        @printf("  greater(%d,0): %+.6f%+.6fi vs ED %+.6f%+.6fi |%.2e\n",
                i, real(gmine), imag(gmine), real(gΔ), imag(gΔ), abs(gmine - gΔ))
    end
    # Z 对照：decoupled Z = Tr[e^{−βH_S}]·Tr[e^{−βh_b}]²（归一化约定无关，只看量级合理性）
    println()
end
println("判决完成")
