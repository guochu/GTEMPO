# reference-solution helpers shared by the model tests. `μ` 是交给模型的化学势
# 约定（H_imp = μ Σ n_i）。ED 构造器（singlemode_ed 等）在顶层 test/util.jl。

using LinearAlgebra: exp, Diagonal

# ED Matsubara Green's function ⟨a(τ) a†⟩ of the full ED Hamiltonian
gτ_ed(H, a, adag, τs, β) = correlation_2op_1τ(H, a, adag, τs, β=β)

"""
    greater_lesser_ed(H, a, adag, H0, ts, β)

ED real-time greater/lesser functions with initial density matrix
`exp(-β H0)` (H0 = the decoupled Hamiltonian for a separable thermal state).
"""
function greater_lesser_ed(H, a, adag, H0, ts, β)
	ρ = exp(-β * H0)
	cache = eigencache(H)
	gt = -im .* correlation_2op_1t(H, a, adag, ρ, ts, cache, reverse=false)
	lt = im .* correlation_2op_1t(H, adag, a, ρ, ts, cache, reverse=true)
	return gt, lt
end

# Kadanoff variant: the τ leg builds the interacting thermal state exp(-β H)
function greater_lesser_ed_mixed(H, a, adag, ts, β)
	ρ = exp(-β * H)
	cache = eigencache(H)
	gt = -im .* correlation_2op_1t(H, a, adag, ρ, ts, cache, reverse=false)
	lt = im .* correlation_2op_1t(H, adag, a, ρ, ts, cache, reverse=true)
	return gt, lt
end

"杂质热初态 e^{−βH_imp}/tr（实时轮廓 initial_state；浴侧真空由 IF 承担）"
thermal_impurity_state(model, β) =
	Matrix{ComplexF64}(fock_thermalstate(model, β).data)

"""
    fermionic_setup(lattice, bath, model; β=nothing, order=1, coupling_norm=nothing)

新 API 装配：逐带 PartialIF（带荷/coupling_norm 随带携带）+ K（可带热初态、
二阶 Weyl 标量）。返回 (K, Is)。
"""
function fermionic_setup(lattice, bath, model; β=nothing, order::Int=1, coupling_norm=nothing)
	order == 2 && coupling_norm === nothing &&
		throw(ArgumentError("order=2 requires coupling_norm"))
	Is = ntuple(b -> hybriddynamics(lattice, bath, PartialIF(); band=b,
	                                coupling_norm=coupling_norm), lattice.bands)
	initial = β === nothing ? nothing : thermal_impurity_state(model, β)
	K = order == 2 ?
		sysdynamics_weyl2order(lattice, model; coupling_norm=coupling_norm, initial_state=initial) :
		sysdynamics(lattice, model; initial_state=initial)
	return K, Is
end
