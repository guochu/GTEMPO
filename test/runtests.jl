using Test, Random
using GTEMPO
import GrassmannMPS as GM
using ImpurityModelBase
# ImpurityModelBase（contour Green 函数容器）也导出 branches，与 GTEMPO 的
# lattice branches 消歧：测试里的 branches 一律指 GTEMPO 的
using GTEMPO: branches

Random.seed!(12354)

# 系统因子 K 构造缺省后端 = :fock（plan/sysdynamics_fock_mpo_plan.md Phase 3）；
# 旧径机制测试已在各自文件内显式钉 :reference。如需全套切回旧径验证：
# GTEMPO.set_sysfactor_backend!(:reference)

# shared utilities: spectra + ED reference builders（API 无关）
include("util.jl")

# Part 1: basic API tests (fast sanity checks of the public interface)
include("api/runtests.jl")

# Part 2: functional tests against exact references (ED / analytic)
include("normalbath/runtests.jl")

# BCS（D7：仅 Z2）与多带（≤2 带）功能测试
include("bcs/runtests.jl")
include("multiorbital/runtests.jl")
