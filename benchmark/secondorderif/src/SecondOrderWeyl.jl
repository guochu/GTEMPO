# ============================================================================
# SecondOrderWeyl — Weyl 二阶 GTEMPO：完整版算法模块
#
# 对齐 U1GTEMPO2（~/Documents/Missile/U1GTEMPO2）的模块结构，但全部建立在
# GTEMPO 已有的 Z2 GMPS 框架之上；对任意 ConstImpurityHamiltonian（任意带数）
# 与一般有界谱浴成立——单元级实时核对 Dirac-δ 谱走闭式精确路径，一般谱走
# QuAPI 单元级积分的谱路径（kernels.jl）：
#
#   U1GTEMPO2 模块       本目录对应          说明
#   ------------------   -----------------   ------------------------------
#   contours.jl          contours.jl         N 个全宽源单元的轮廓（无端点）
#   kernels.jl           kernels.jl          Weyl N×N 浴核 + coupling_norm
#   exponentials.jl      （并入 kernels）    单元积分闭式（Dirac-δ）/ 谱路径（一般浴）
#   tensors.jl           tensors.jl          单元核 → Z2 site tensors
#   system.jl            system.jl           Weyl/一阶系统因子 K（folded 链）
#   influence.jl         influence.jl        IF：完全复用 GTEMPO hybriddynamics
#   contraction.jl       contraction.jl      积分：完全复用 GTEMPO integrate
#   solver.jl            solver.jl           Problem/Solution/Solve
#   measurement.jl       measurement.jl      occupation / greater / lesser
#   prepared_system.jl   （并入 system）     基准规模无需环境缓存
#   scaling.jl           （缩放折叠进张量）  GTEMPO GMPS 自带 scaling 字段
#   persistence.jl       （基准不需要）
#
# 与 U1GTEMPO2 的核心差异：U1GTEMPO2 的收缩实现 Gaussian 测度，其 Weyl 局部
# 权重做了测度适配（aā 系数 = I/2+n）；GTEMPO 的 integrate 是裸 Berezin 配对，
# 相干态高斯权重在一阶管线里放在传播子的 overlap 因子中（bulkconnection）。
# Weyl K 是自洽的 automaton 链（Fock 键穿线、无 overlap、无边界闭合），因此
# 这里直接用推导的裸测度权重（aā 系数 = I），见 system.jl 的注释。
# ============================================================================
module SecondOrderWeyl

using GTEMPO
using GTEMPO: GrassmannContourIndex, GrassmannMPS, GrassmannLattice,
              RealGrassmannLattice, ImagGrassmannLattice,
              z2space, vacuumstate, hybriddynamics, integrate,
              truncdimcutoff, TruncationScheme, swapband,
              ExactIntegrate, fock_propagator, unset_svectors!
using Z2Tensors
using Z2Tensors: Z2Space, Z2Irrep, TensorMap, ←
using LinearAlgebra: eigen, Diagonal, tr, norm, svd, dot, I, diagm
using QuadGK: quadgk

const _PH = Z2Space(0 => 1, 1 => 1)   # 单个 Grassmann 变量的物理空间

include("contours.jl")
include("kernels.jl")
include("tensors.jl")
include("system.jl")
include("influence.jl")
include("contraction.jl")
include("solver.jl")
include("measurement.jl")
include("diagnostics.jl")

# —— 轮廓（对齐 U1GTEMPO2 contours.jl）
export AbstractWeylContour, WeylRealTime, WeylImagTime, WeylMixedTime
export weyl_lattice, ncells, stepλ, cell_interval, boundary_times, cell_indices
# —— 浴核（对齐 kernels.jl / exponentials.jl）
export HybridizationKernel, hybridization, coupling_norm
export realtime_kernel_trace, realtime_kernel_quapi, α_weight
# —— 系统因子（对齐 system.jl）
export local_weight, system_factor, WeylInsertion
# —— IF / 收缩 / 求解 / 观测量（对齐 influence / contraction / solver / measurement）
export influence_functional
export partition_function, expectation
export Problem, Solution, solve
export occupation, greater, lesser, raw_correlation, matsubara
# —— 诊断
export fit_order, convergence_report, maxbond
# —— 内部（校准/调试用）
export CellCore, split_core, _fock_annihilators, KernelCorrView

end # module
