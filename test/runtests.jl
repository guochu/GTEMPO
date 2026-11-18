using Test, Random
using GTEMPO
import GrassmannMPS as GM
using ImpurityModelBase
# ImpurityModelBase（contour Green 函数容器）也导出 branches，与 GTEMPO 的
# lattice branches 消歧：测试里的 branches 一律指 GTEMPO 的
using GTEMPO: branches

Random.seed!(12354)

# 系统因子 K 构造 = fock_sysdynamics + fockmpo2gmps 唯一路径（旧径
# SystemFactorPreparation/_system_factor_reference 已随 V2 删除）。

# shared utilities: spectra + ED reference builders（API 无关）
include("util.jl")

# Part 1: basic API tests (fast sanity checks of the public interface)
include("api/runtests.jl")

# Part 2: functional tests against exact references (ED / analytic)
include("normalbath/runtests.jl")

# BCS（D7：仅 Z2）与多带（≤2 带）功能测试
include("bcs/runtests.jl")
include("multiorbital/runtests.jl")
