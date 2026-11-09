# =============================================================================
# toulouse/realtime.jl —— Toulouse 模型（U=0, ε_d=0.7）基于严格 IF 的 1/2/4 阶
# Trotter 实时演化 error scaling（紧凑消元路线；替换旧的单扇区 gt/lt 诊断脚本
# ——旧脚本受单扇区闭合约定差异限制，order 4 不收敛，已删除）。
#
# 判据：
#   ① 机器级预检：紧凑管线 n(T) == Fock 离散 Trotter 链（同 schedule、同热态 ρ₀）；
#   ② error scaling：|n_N − n_exact| vs δt = T/N（n_exact = 4 维全 H 精确对角化）；
#   ③ 真阶判定：order 1/2 → p≈2（对角 H_imp 下严格核精确重求和 [A,B]δt²），
#     order 4 → p≈4（每减半 δt 误差 ÷16）。
#
# 输出：本文件夹 realtime.csv / realtime.png。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/toulouse/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "commons", "scaling_compact_models.jl"))
run_realtime_scaling("Toulouse")
