# impuritymodels：杂质模型层入口
#
# 保留（零改动）：predefined 模型（AndersonIM/ToulouseIM/IRLM/KanamoriIM）与
# generalimpurity 的 Fock 空间纯数值接口（fockmatrix/fock_propagator/
# fock_thermalstate——cell automaton K 构造的数据源，src/sysdynamics.jl 消费）。
# 旧装配路径（_propagator_sparsegmps 窗口 mult!、sysdynamics_fast、bare 族、
# 8 个 *_deprecated、termdynamics、fermionoperator）已随 D2 裁决删除；
# K 构造统一走 cell automaton（Phase 2 接线）。

"""
	AbstractImpurityHamiltonian

Supertype of all impurity Hamiltonians. Interface: `fock_propagator`
(both subtypes) and `fock_thermalstate` (needed to set the thermal
initial state on real-time lattices); the number of impurity bands is
reported by `num_bands`.
"""
abstract type AbstractImpurityHamiltonian end
"""
	ConstImpurityHamiltonian <: AbstractImpurityHamiltonian

Constant (time-independent) impurity Hamiltonian; quench Hamiltonians
belong here. Interface: `fock_propagator(model, branch, dt)`.
"""
abstract type ConstImpurityHamiltonian <: AbstractImpurityHamiltonian end
"""
	AbstractTdImpurityHamiltonian <: AbstractImpurityHamiltonian

Time-dependent impurity Hamiltonian. Interface: `fock_propagator(model,
branch, dt, t)` on the real-time branches (`:+`/`:-`), with `t` the left
endpoint of the physical time interval of the step.
"""
abstract type AbstractTdImpurityHamiltonian <: AbstractImpurityHamiltonian end

# interface for AbstractImpurityHamiltonian
fock_propagator(model::ConstImpurityHamiltonian, branch::Symbol, dt::Real) =
    error("fock_propagator not implemented for model $(typeof(model))")
fock_thermalstate(h::AbstractImpurityHamiltonian, β::Real) =
    error("fock_thermalstate not implemented for model $(typeof(h))")
fock_propagator(model::AbstractTdImpurityHamiltonian, branch::Symbol, dt::Real, t::Real) =
    error("fock_propagator not implemented for model $(typeof(model))")
# constant models: the propagator on the real-time branches is time-independent,
# the branch time t is ignored
fock_propagator(model::ConstImpurityHamiltonian, branch::Symbol, dt::Real, t::Real) =
    fock_propagator(model, branch, dt)

"""
	num_bands(h) -> Int

Number of impurity bands of the model. Defaults to the `bands` field of
the model; models without such a field define their own method.
"""
num_bands(h::AbstractImpurityHamiltonian) = h.bands

include("predefined/siam.jl")
include("predefined/irlm.jl")
include("predefined/skmodel.jl")

include("generalimpurity/generalimpurity.jl")

# 公开皮肤：杂质模型类型与接口（旧 export 名保持）
export AbstractImpurityHamiltonian, ConstImpurityHamiltonian, AbstractTdImpurityHamiltonian,
       ImpurityHamiltonian, QuenchedImpurityHamiltonian, TdImpurityHamiltonian, TdImpurityOp,
       AndersonIM, ToulouseIM, IRLM, KanamoriIM
export num_bands, fock_propagator, fock_thermalstate, fockmatrix, fockstate
