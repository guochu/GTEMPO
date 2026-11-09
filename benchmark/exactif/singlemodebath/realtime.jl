# =============================================================================
# singlemodebath/realtime.jl —— 单模 bath Anderson 模型（μ=0.7, U=1.5, bands=2）
# 基于严格 IF 的 1/2/4 阶 Trotter 实时演化 error scaling（紧凑消元路线；替换
# 旧的单扇区 gt/lt 诊断脚本——旧脚本受单扇区闭合约定差异限制，order 4 不收敛，
# 已删除）。
#
# 判据：
#   ① 机器级预检：紧凑管线 n(T) == Fock 离散 Trotter 链（同 schedule、同热态
#     ρ₀；16 维 Fock 直算，含 U 关联）；
#   ② error scaling：|n_N − n_exact| vs δt = T/N（n_exact = 16 维全 H 精确对角化）；
#   ③ 真阶判定：order 1/2 → p≈2、order 4 → p≈4（每减半 δt 误差 ÷16；
#     粗端 δt=0.4 处 order 4 仅略优于 order 1 属 Suzuki 负 cell 预渐近，预期行为）。
#
# 输出：本文件夹 realtime.csv / realtime.png。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/singlemodebath/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "commons", "scaling_compact_models.jl"))
run_realtime_scaling("Anderson")
