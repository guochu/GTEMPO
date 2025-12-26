# arXiv:1407.1832v1 "Time-evolving a matrix product state with long-ranged interactions"
abstract type TimeEvoMPOAlgorithm <: MPSAlgorithm end
abstract type FirstOrderStepper <: TimeEvoMPOAlgorithm end
abstract type SecondOrderStepper <: TimeEvoMPOAlgorithm end

@with_kw struct WI <: FirstOrderStepper
	tol::Float64 = Defaults.tol
	maxiter::Int = Defaults.maxiter
end

@with_kw struct WII <: FirstOrderStepper
	tol::Float64 = Defaults.tol
	maxiter::Int = Defaults.maxiter
end

struct ComplexStepper{F<:FirstOrderStepper} <: SecondOrderStepper
	stepper::F
end

function get_A(x::SchurMPOTensor)
	s1, s2 = size(x)
	return x.Os[2:s1-1, 2:s2-1]
end
function get_B(x::SchurMPOTensor)
	s1, s2 = size(x)
	return x.Os[2:s1-1, end]
end
function get_C(x::SchurMPOTensor)
	s1, s2 = size(x)
	return x.Os[1, 2:s2-1]
end
get_D(x::SchurMPOTensor) = x[1, end]

function _SiteW_impl(WA, WB, WC, WD)
	s1, s2 = size(WA)
	r = Matrix{Any}(undef, s1+1, s2+1)
	r[1, 1] = WD
	for l in 2:s2+1
		r[1, l] = WC[l-1]
	end
	for l in 2:s1+1
		r[l, 1] = WB[l-1]
	end
	r[2:end, 2:end] .= WA
	return SparseMPOTensor(r)
end

function _sqrt2(dt::Complex) 
	r = sqrt(dt)
	return r, r
end

function _sqrt2(dt::Real)
	if dt >= zero(dt)
	 	r = sqrt(dt)
	 	return r, r
	 else
	 	r = sqrt(-dt)
	 	return r, -r
	end 
end

function timeevompo(m::SchurMPOTensor, dt::Number, alg::WI)
	WA = get_A(m)
	δ₁, δ₂ = _sqrt2(dt)
	WB = get_B(m) .* δ₁
	WC = get_C(m) .* δ₂
	D = get_D(m)
	WD = isomorphism(storagetype(D), codomain(D), domain(D)) + dt * D
	return _SiteW_impl(WA, WB, WC, WD)
end


function timeevompo(m::SchurMPOTensor, dt::Number, alg::WII)
	s1, s2 = size(m)
	_S = promote_type(scalartype(m), scalartype(dt))
	T = mpotensortype(spacetype(m), _S)
	WA = Matrix{T}(undef, s1-2, s2-2)
	WB = Vector{T}(undef, s1-2)
	WC = Vector{T}(undef, s2-2)
	D = convert(T, get_D(m)) 
	WD = D

	δ = dt
	δ₁, δ₂ = _sqrt2(dt)
	for j in 2:s1-1, k in 2:s2-1
		init_1 = isometry(storagetype(D), codomain(D), domain(D))
		init = [init_1, convert(T, zero(m[1, k])), convert(T, zero(m[j, end])), convert(T, zero(m[j, k]))]

		(y, convhist) = exponentiate(1.0, init, Arnoldi(; tol=alg.tol, maxiter=alg.maxiter)) do x
			out = similar(x)

			@tensor out[1][-1 -2; -3 -4] := δ * x[1][-1 1; -3 -4] * m[1, end][2 -2; 2 1] 

			@tensor out[2][-1 -2; -3 -4] := δ * x[2][-1 1; -3 -4] * m[1, end][2 -2; 2 1] 

			@tensor out[2][-1 -2; -3 -4] += δ₂ * x[1][4 3; 4 -4] * m[1, k][-1 -2; -3 3] 

			@tensor out[3][-1 -2; -3 -4] := δ * x[3][-1 1; -3 -4] * m[1, end][2 -2; 2 1] 

			@tensor out[3][-1 -2; -3 -4] += δ₁ * x[1][4 3; 4 -4] * m[j, end][-1 -2; -3 3]

			@tensor out[4][-1 -2; -3 -4] := δ * x[4][-1 1; -3 -4] * m[1, end][2 -2; 2 1] 

			@tensor out[4][-1 -2; -3 -4] += x[1][4 3; 4 -4] * m[j, k][-1 -2; -3 3]  

			@tensor out[4][-1 -2; -3 -4] += δ₁ * x[2][4 3; -3 -4] * m[j, end][-1 -2; 4 3]  

			@tensor out[4][-1 -2; -3 -4] += δ₂ * x[3][-1 4; 3 -4] * m[1, k][3 -2; -3 4] 

			return out
		end
		convhist.converged == 0 && @warn "failed to exponentiate $(convhist.normres)"

		WA[j - 1, k - 1] = y[4]
		WB[j - 1] = y[3]
		WC[k - 1] = y[2]
		WD = y[1]

	end
	return _SiteW_impl(WA, WB, WC, WD)
end

# function _W_impl(h::Vector)
# 	h2 = copy(h)
# 	# h2[1] = SparseMPOTensor(h[1].Os[1:1, :], h[1].leftspaces[1:1], h[1].rightspaces, h[1].pspace)
# 	# h2[end] = SparseMPOTensor(h[end].Os[:, 1:1], h[end].leftspaces, h[end].rightspaces[1:1], h[end].pspace)
# 	h2[1] = h[1][1,:]
# 	h2[end] = h[end][:,1]
# 	return MPOHamiltonian(h2)
# end

timeevompo(m::MPOHamiltonian{<:SchurMPOTensor}, dt::Number, alg::FirstOrderStepper) = MPOHamiltonian([timeevompo(mj, dt, alg) for mj in m.data])
function timeevompo(h::Union{SchurMPOTensor, MPOHamiltonian{<:SchurMPOTensor}}, dt::Number, alg::ComplexStepper)
	dt1, dt2 = complex_stepper(dt)
	return timeevompo(h, dt1, alg.stepper), timeevompo(h, dt2, alg.stepper)
end
timeevompo(m::MPOHamiltonian{<:SchurMPOTensor}, dt::Number; alg::TimeEvoMPOAlgorithm = WII()) = timeevompo(m, dt, alg)

"""
	complex_stepper(dt::Number)
Retun dt₁, dt₂
If U = exp(H*dt) is a first order stepper,
then U₁ = exp(H*dt₁), U₂ = exp(H*dt₂), and
U = U₁U₂ will be a second order stepper
"""
complex_stepper(dt::Number) = (1-im) * dt/2, (1+im) * dt/2
