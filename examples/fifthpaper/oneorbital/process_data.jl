# push!(LOAD_PATH, "../../../src")

using ImpurityModelBase
using JSON, QuadGK


parse_complex_number(a) = a["re"] + im * a["im"]

function read_itempo(t; U=1, μ=U/2, order=6, chi=1024)
	U = convert(Float64, U)
	μ = convert(Float64, μ)
	dt = 0.05
	k = 5
	prony = 5
	N = round(Int, t / dt)
	file_name = "result/SIAM_onebath_beta40.0_U$(U)_mu$(μ)_dt$(dt)_k$(k)_trunc$(order)_prony$(prony)_N$(N)_chi$(chi).json"
	data = JSON.parsefile(file_name)

	gt = [im*parse_complex_number(item) for item in data["gt"]]
	lt = [im*parse_complex_number(item) for item in data["lt"]]
	# println(gt[1:10])
	# println(lt[1:10])

	# error("stop here")
	retarded = gt + lt

	return -retarded 
end

# function compute_Gτ(retarded_t, t, lb, ub)
# 	β = 40.
# 	dt = 0.05
# 	r = LinearPrediction(retarded_t, stepsize=dt)
# 	# t = 100.
# 	n = round(Int, t/dt)
# 	retarded = [r[i] for i in 1:n]

# 	# println(retarded[end-10:end])

# 	# m = 10
# 	# dt = dt / m
# 	# # retarded, dt = scaling_down(retarded, dt, m)
# 	# retarded = r.(0:dt:200)

# 	Gw = FourierTransform(retarded, δt=dt, δ=1.0e-8)
# 	Aw(ϵ) = begin 
# 		x = -imag(Gw(ϵ))/π
# 		return (x >= 0) ? x : zero(x)
# 		# return x
# 	end

# 	# println(quadgk(ϵ -> Aw(ϵ), -1, 2.))

# 	# return Aw

# 	# lb = -2
# 	# ub = 2
# 	nm = quadgk(ϵ -> Aw(ϵ), lb, ub)[1]
# 	println("norm is ", nm)

# 	f(ϵ, τ) = -exp(-ϵ*τ) * Aw(ϵ)/((1+exp(-β*ϵ)) * nm)
# 	G(τ) = quadgk(ϵ -> f(ϵ, τ), lb, ub)[1]

# 	δτ = 0.1
# 	τs = collect(0:δτ:β)
# 	ws = collect(lb:0.05:ub)
# 	return ws, Aw.(ws), τs, G.(τs)
# end

function compute_Gτ(retarded_t, lb, ub)
	println("final retarded_t value ", retarded_t[end])
	β = 40.
	δt = 0.001
	retarded = linear_predict(retarded_t, 0.05, δt=δt)

	println("final retarded value ", retarded[end])
	dw = 1.0e-4
	ws = collect(lb:dw:ub)
	Aw = Gt_to_Aw(retarded, δt, lb=lb, ub=ub, dw=dw, verbosity=2)

	δτ = 0.1
	τs = collect(0:δτ:β)
	Gtau = Aw_to_Gτ(Aw, β=β, lb=lb, ub=ub, dw=dw, δτ=δτ)

	# ws2 = collect(lb:0.05:ub)
	# Aw2 = Gt_to_Aw(retarded, δt, lb=lb, ub=ub, dw=0.05)
	return ws, Aw, τs, Gtau
end

function main(t; U=1, μ=U/2, lb=-2, ub=2, order=6, chi=1024)
	lb = convert(Float64, lb)
	ub = convert(Float64, ub)
	t = convert(Float64, t)
	U = convert(Float64, U)
	μ = convert(Float64, μ)

	dt = 0.05
	k = 5
	prony = 5
	N = round(Int, t / dt)


	retarded = read_itempo(t, U=U, μ=μ, order=order, chi=chi)
	ws, Aw, τs, Gτ = compute_Gτ(retarded, lb, ub)

	file_name = "result/SIAM_onebath_beta40.0_t$(t)_U$(U)_mu$(μ)_lb$(lb)_ub$(ub)_dt$(dt)_k$(k)_trunc$(order)_prony$(prony)_chi$(chi).json"

	results = Dict("ts"=>τs, "gf"=>Gτ, "ws"=>ws, "Aw"=>Aw)

	println("save results to ", file_name)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return results
end

function main_all_U(t = 20.; lb=-2, ub=2, order=6, chi=1024)
	for U in 0.1:0.1:1.
		main(t; U=U, lb=lb, ub=ub, order=order, chi=chi)
	end
end

