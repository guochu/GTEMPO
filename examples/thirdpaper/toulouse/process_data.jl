push!(LOAD_PATH, "../../../src")

using GTEMPO, ImpurityModelBase.Utilities
using JSON, QuadGK
# using Interpolations


parse_complex_number(a) = a["re"] + im * a["im"]

# function compute_GW(data)
# 	# println(data)
# 	greater = [parse_complex_number(item) for item in data[1]]
# 	lesser = [parse_complex_number(item) for item in data[2]]

# 	obs = greater - lesser

# 	# dt = 0.05
	
# 	# r = LinearExtrapolation(obs, stepsize=dt, nfit=150)

# 	# ts = 0:0.01:50
# 	# gf = [r(t) for t in ts]

# 	# f = FourierTransform(gf, δt=0.01)

# 	f = FourierTransform(obs, δt=0.01)

# 	ws = -10:0.05:10

# 	return ws, -imag(f.(ws))
# end

function read_data(t; ϵ_d=1., δt=0.05, order=7, chi=1024)
	t = convert(Float64, t)
	ϵ_d = convert(Float64, ϵ_d)
	β = 40.

	file_name = "result/thouless_tempo_real_beta$(β)_t$(t)_mu$(ϵ_d)_dt$(δt)_order$(order)_chi$(chi).json"
	data = JSON.parsefile(file_name)

	ns = data["ns"]
	gt = [-im*parse_complex_number(item) for item in data["gt"]]
	lt = [im*parse_complex_number(item) for item in data["lt"]]
	retarded = gt + lt

	# dt = 0.05
	# r = LinearExtrapolation(retarded, stepsize=dt)

	# f = FourierTransform(retarded, δt=0.05)
	return ns, retarded 
end

function compute_Gτ(retarded_t, t, lb, ub)
	β = 40.
	dt = 0.05
	r = LinearPrediction(retarded_t, stepsize=dt)
	# t = 100.
	n = round(Int, t/dt)
	retarded = [r[i] for i in 1:n]

	# println(retarded[end-10:end])

	# m = 10
	# dt = dt / m
	# # retarded, dt = scaling_down(retarded, dt, m)
	# retarded = r.(0:dt:200)

	Gw = FourierTransform(retarded, δt=dt, δ=1.0e-8)
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
	return ws, Aw.(ws), τs, G.(τs)
end

function main(t; t2=100, ϵ_d=1., lb=-2, ub=2, order=7, chi=1024)
	t2 = convert(Float64, t2)
	lb = convert(Float64, lb)
	ub = convert(Float64, ub)
	t = convert(Float64, t)
	ϵ_d = convert(Float64, ϵ_d)
	β = 40.
	δt = 0.05
	ns, retarded = read_data(t, ϵ_d=ϵ_d, δt=δt, order=order, chi=chi)
	ws, Aw, τs, Gτ = compute_Gτ(retarded, t2, lb, ub)

	file_name = "result/thouless_tempo_real_beta$(β)_t$(t)_mu$(ϵ_d)_tep$(t2)_lb$(lb)_ub$(ub)_dt0.05_order$(order)_chi$(chi).json"

	results = Dict("ns" => ns, "ts"=>τs, "gf"=>Gτ, "ws"=>ws, "Aw"=>Aw)

	println("save results to ", file_name)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return results
end
