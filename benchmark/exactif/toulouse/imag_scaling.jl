# =============================================================================
# toulouse/imag_scaling.jl —— Toulouse 模型虚轴 Trotter error scaling
# （orders 1/2/4 × N 扫描；β=1, εd=0.7, WB=1, AB=0.5）
#
# 误差 = |线性调度离散链锚 − 连续极限 ED|。线性调度保留端点半门
# （环合并会破坏 order-4 的 Suzuki 条件——见 imag_theory_notes.md）。
# GMPS 单 MPS 管线对离散链逐位复现（机器精度验证），故此表即管线误差。
#
# 输出：本文件夹 imag_scaling.csv；图：../commons/plot_imag_scaling.py 生成
# 本文件夹 imag_scaling.png。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/toulouse/imag_scaling.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "NativeExact.jl"))
include(joinpath(@__DIR__, "..", "commons", "imag_commons.jl"))
using .ImagCommons
using .NativeExact
using Printf

const MODEL = ModelParams(tag = "toulouse", WB = 1.0, AB = 0.5,
                          BETA = 1.0, EPSD = 0.7)
const ORDERS = [1, 2, 4]
const NS = [1, 2, 4, 8, 16, 32, 64]

ops = build_model(MODEL)
schedfn = (order, N, dt) -> NativeExact.native_schedule(order, N, dt)
zc, nc = scaling_table(ops, MODEL, schedfn,
                       joinpath(@__DIR__, "imag_scaling.csv"), ORDERS, NS)
print_header(MODEL, zc, nc)
@printf("已写 imag_scaling.csv（orders %s × N %s）\n",
        join(ORDERS, ","), join(NS, ","))
