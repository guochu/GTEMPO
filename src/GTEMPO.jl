module GTEMPO

# =====================================================================
# GTEMPO：GrassmannMPS 后端（换轨版，plan/gtempo_grassmannmps_backend_plan.md）
#
# 旧栈（Z2Tensors + src/grassmanntensor + src/grassmannmps）整体退役；
# 张量引擎 = 独立包 GrassmannMPS（TensorKit sector 栈）。业务 API 皮肤
# （模型 / lattice / IF / observables 公开函数名）保持不变。
#
# 功能调整（§0）：
#   #1 GrassmannOrdering 家族与 Chrono ordering 类型层删除（排布固定 Chrono，
#       lattice 无 ordering 输入——2026-11-16）；
#   #2 electronphonon 暂停（整目录移除，恢复路径见 changes.md）；
#   #3 sysdynamics_weyl2order 新增（port/sysdynamics.jl，K 链局部 Weyl 二阶）；
#   #4 对称性用户手选 :z2/:z2z2/:u1/:u1u1（port/bridge.jl，D1 裁决）；
#      BCS 仅 :z2（D7 裁决）；CUDA ext 直接删除（用户裁决）。
# =====================================================================

using Base: @boundscheck, @propagate_inbounds
using Reexport
@reexport using ImpurityModelBase, QuAPI
import QuAPI: branch
using LinearAlgebra: LinearAlgebra, eigen, Diagonal, I, tr, norm, diag
import TensorKit as TK
import MatrixAlgebraKit as MAK
using TensorKit: ⊗, ⊠, ←, Vect, Z2Irrep, U1Irrep, Sector,
                 sectors, sectortype, dim, fusiontrees, SpaceMismatch,
                 SectorVector, SectorDict
# GrassmannMPS 后端（R7：选择性 import，勿整包 using——norm/dot/lmul! 与
# LinearAlgebra 的重载关系按 GM 模块约定处理；R8：私有原语只 import 不扩方法）
import GrassmannMPS as GM
import GrassmannMPS: vacuumstate, expectation,   # 扩展目标（lattice / solution 方法）
    coefficient_scaled                           # 桥层扩展目标（同函数新方法）
using GrassmannMPS: GMPS, SparseGMPS, ScaledScalar, ExpectationCache,
    ChargeCharge, Z2Z2, grassmannspace, grassmann_basis, wedge_vertex,
    canonical_grassmann_spaces,
    environments, Zvalue, value, truncrelerr, truncdimcutoff, truncate!,
    leftorth!, rightorth!, swap!, reorder!, randomgmps,
    mult, mult!, apply!, DMRG2, SVDCompression,
    canonicalize!, iscanonical, entanglement_entropy, entanglement_spectrum,
    schmidt_values, bonddims, bondspaces, physpace, physpaces,
    logscale, setlogscale!, todense, coefficient_scaled, coefficient, NoTruncation

# ---------------------------------------------------------------------
# 后端适配层（GrassmannMPS 桥 + Chrono lattice + 缩放）+ 相关函数层
# ---------------------------------------------------------------------
include("bridge.jl")
include("lattices.jl")

# correlation functions（QuAPI 薄壳，零逻辑改动）
include("correlationfunction.jl")

# 杂质模型（fock_propagator/fock_thermalstate 纯数值接口）
include("impuritymodels/impuritymodels.jl")

# K 构造（cell automaton；sysdynamics / sysdynamics_weyl2order）
include("sysdynamics.jl")

# IF 构造（PartialIF / ExactTTI 首批，D3）
include("influence.jl")

# 积分与观测量（GM 引擎接线 + inserted-K 路线）
include("measurement.jl")

# ---------------------------------------------------------------------
# export 皮肤
# ---------------------------------------------------------------------

# auxiliary / 截断
export truncdimcutoff, truncrelerr, NoTruncation

# GrassmannMPS 皮肤（旧 export 名 → GM 对应物；plan §Phase 1 对表）
export scalartype, @grassmann
export GrassmannMPS, SparseGMPS, AbstractGMPS, togmps, randomgmps
export mult, mult!, apply!, DMRG2, SVDCompression
export leftorth!, rightorth!, canonicalize!, iscanonical
export swap!, reorder!
export entanglement_entropy, entanglement_spectrum, schmidt_values
export bonddims, bondspaces, physpace, physpaces, logscale, setlogscale!

# lattice（Chrono 排布 + 对称性选择器）
export AbstractGrassmannLattice, AbstractImagGrassmannLattice, AbstractRealGrassmannLattice,
       AbstractMixedGrassmannLattice
export ImagGrassmannLattice, RealGrassmannLattice, MixedGrassmannLattice, GrassmannLattice
export branches, indexmappings, vacuumstate
export symmetry, site_spaces

# 对称性 / 桥层（构造 + 引擎接线）
export Truncation
export pairize, monomial, truncate!
export charge_dimensions, check_charge
export coefficient, coefficient_scaled, integrate, integrate_scaled
export EnvironmentCache
export ScaledScalar, value

# correlation functions
export correlationfunction

# K 构造
export sysdynamics, sysdynamics_weyl2order, OperatorInsertion

# IF 构造
export InfluenceFunctionalAlgorithm, PartialIF, ExactTTI, InfluenceFunctional
export bcs_hybriddynamics, bcs_partial_if
export symmetry_spec, var_charge, cell_charges, fock_charges
export HybridizationKernel, hybridization, influence_functional, hybriddynamics

# 观测量（Phase 5b：inserted-K 边界观测 + 算符缓存原始关联两轨，D5 裁决见 measurement.jl）
export occupation, greater, lesser, matsubara, default_branch, nn
export lattice_monomial, environments_cache, cached_Z, cached_gf, cached_correlation

end # module
