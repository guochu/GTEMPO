push!(LOAD_PATH, "../../../src")

using ImpurityModelBase
using JSON, QuadGK
using LinearAlgebra


parse_complex_number(a) = a["re"] + im * a["im"]

function read_data(t, t₀; U=1, μ=U/2, order=6, chi=1024)
	t = convert(Float64, t)
	t₀ = convert(Float64, t₀)
	U = convert(Float64, U)
	μ = convert(Float64, μ)
	file_name = "result/anderson_tempo1_beta40.0_t$(t)_$(t₀)_U$(U)_mu$(μ)_dt0.05_order$(order)_chi$(chi).json"
	data = JSON.parsefile(file_name)

	ns = data["ns"]
	gt = [-im*parse_complex_number(item) for item in data["gt"]]
	lt = [im*parse_complex_number(item) for item in data["lt"]]
	retarded = gt + lt

	return ns, retarded 
end


function compute_Gτ(retarded_t, t, lb, ub)
	println("final retarded_t value ", retarded_t[end])
	β = 40.
	dt = 0.05
	# r = LinearPrediction(retarded_t, stepsize=dt, nfit=div(length(retarded_t), 2))
	r = LinearPrediction(retarded_t, stepsize=dt)
	# t = 100.
	# n = round(Int, t/dt)
	# retarded = [r[i] for i in 1:n]

	dt2 = 0.01
	n = round(Int, t/dt2)
	retarded = [r(i*dt2) for i in 1:n]


	println("final retarded value ", retarded[end])
	# retarded = im * imag(retarded)

	# println(retarded[end-10:end])

	# m = 10
	# dt = dt / m
	# # retarded, dt = scaling_down(retarded, dt, m)
	# retarded = r.(0:dt:200)

	Gw = FourierTransform(retarded, δt=dt2, δ=1.0e-8)
	Aw(ϵ) = begin 
		x = -imag(Gw(ϵ))/π
		return (x >= 0) ? x : zero(x)
		# return x
	end

	# println(quadgk(ϵ -> Aw(ϵ), -1, 2.))

	# return Aw

	# lb = -2
	# ub = 2
	nm = quadgk(ϵ -> Aw(ϵ), lb, ub)[1]
	println("norm is ", nm)

	f(ϵ, τ) = -exp(-ϵ*τ) * Aw(ϵ)/((1+exp(-β*ϵ)) * nm)
	G(τ) = quadgk(ϵ -> f(ϵ, τ), lb, ub)[1]

	δτ = 0.1
	τs = collect(0:δτ:β)
	ws = collect(lb:0.05:ub)
	return ws, Aw.(ws) ./  nm, τs, G.(τs)
end

function compute_Gτ_2(retarded_t, t, lb, ub)
	println("final retarded_t value ", retarded_t[end])
	β = 40.
	δt = 0.0001
	retarded = linear_predict(retarded_t, 0.05, δt=δt)

	println("final retarded value ", retarded[end])
	dw = 1.0e-4
	ws = collect(lb:dw:ub)
	Aw = Gt_to_Aw(retarded, δt, lb=lb, ub=ub, dw=dw, verbosity=2)

	δτ = 0.1
	τs = collect(0:δτ:β)
	Gtau = Aw_to_Gτ(Aw, β=β, lb=lb, ub=ub, dw=dw, δτ=δτ)

	ws2 = collect(lb:0.05:ub)
	Aw2 = Gt_to_Aw(retarded, δt, lb=lb, ub=ub, dw=0.05)
	return ws2, Aw2, τs, Gtau
end

function main(t, t₀=t/2; t2=100, U=1, μ=U/2, lb=-2, ub=2, order=6, chi=1024)
	t2 = convert(Float64, t2)
	t₀ = convert(Float64, t₀)
	lb = convert(Float64, lb)
	ub = convert(Float64, ub)
	t = convert(Float64, t)
	U = convert(Float64, U)
	μ = convert(Float64, μ)
	ns, retarded = read_data(t, t₀, U=U, μ=μ, order=order, chi=chi)
	@time ws, Aw, τs, Gτ = compute_Gτ(retarded, t2, lb, ub)

	@time ws2, Aw2, τs2, Gτ2 = compute_Gτ_2(retarded, t2, lb, ub)

	# return ws, Aw, τs, Gτ, ws2, Aw2, τs2, Gτ2
	return norm(Aw - Aw2) / norm(Aw), norm(Gτ - Gτ2) / norm(Gτ)
end

