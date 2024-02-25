push!(LOAD_PATH, "../../../src")
push!(LOAD_PATH, "../../../../FermionicTCMPS/src")

using GTEMPO, FermionicTCMPS.Utilities
using JSON, QuadGK
# using Interpolations


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

function read_data_2(t, t₀; U=1, μ=U/2, order=6, chi=1024)
	t = convert(Float64, t)
	t₀ = convert(Float64, t₀)
	U = convert(Float64, U)
	μ = convert(Float64, μ)
	file_name = "result/anderson_tempo1_beta40.0_t$(t)_$(t₀)_U$(U)_mu$(μ)_dt0.05_order$(order)_chi$(chi)_2.json"
	data = JSON.parsefile(file_name)

	ns = data["ns"]
	gt = [-im*parse_complex_number(item) for item in data["gt"]]
	lt = [im*parse_complex_number(item) for item in data["lt"]]
	retarded = gt + lt

	return ns, retarded 
end

# function scaling_down(retarded, dt, m)
# 	dt′ = dt / m
# 	# ts = [i*dt for i in 0:length(retarded)-1]
# 	ts = 0:dt:(length(retarded)-1)*dt
# 	# f1 = cubic_spline_interpolation(ts, real(retarded))
# 	# f2 = cubic_spline_interpolation(ts, imag(retarded))
# 	f1 = linear_interpolation(ts, real(retarded))
# 	f2 = linear_interpolation(ts, imag(retarded))

# 	retarded1 = f1.(0:dt′:ts[end])
# 	retarded2 = f2.(0:dt′:ts[end])
# 	retarded′ = retarded1 + im .* retarded2
# 	return retarded′, dt′
# end


function compute_Gτ(retarded_t, t, lb, ub)
	β = 40.
	dt = 0.05
	r = LinearExtrapolation(retarded_t, stepsize=dt)
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

function main(t, t₀=t/2; t2=100, U=1, μ=U/2, lb=-2, ub=2, order=6, chi=1024)
	t2 = convert(Float64, t2)
	t₀ = convert(Float64, t₀)
	lb = convert(Float64, lb)
	ub = convert(Float64, ub)
	t = convert(Float64, t)
	U = convert(Float64, U)
	μ = convert(Float64, μ)
	ns, retarded = read_data(t, t₀, U=U, μ=μ, order=order, chi=chi)
	ws, Aw, τs, Gτ = compute_Gτ(retarded, t2, lb, ub)

	file_name = "result/anderson_tempo1_beta40.0_t$(t)_U$(U)_mu$(μ)_tep$(t2)_lb$(lb)_ub$(ub)_dt0.05_order$(order)_chi$(chi).json"

	results = Dict("ns" => ns, "ts"=>τs, "gf"=>Gτ, "ws"=>ws, "Aw"=>Aw)

	println("save results to ", file_name)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return results
end

function main_all_U(t₀, t = t₀ + 20.; t2=100, lb=-2, ub=2, order=6, chi=1024)
	for U in [0.1, 0.5, 1.]
		main(t, t₀; t2=t2, U=U, lb=lb, ub=ub, order=order, chi=chi)
	end
end

function main2(t, t₀=t/2; t2=100, U=1, μ=U/2, lb=-2, ub=2, order=6, chi=1024)
	t2 = convert(Float64, t2)
	t₀ = convert(Float64, t₀)
	lb = convert(Float64, lb)
	ub = convert(Float64, ub)
	t = convert(Float64, t)
	U = convert(Float64, U)
	μ = convert(Float64, μ)
	ns, retarded = read_data_2(t, t₀, U=U, μ=μ, order=order, chi=chi)
	ws, Aw, τs, Gτ = compute_Gτ(retarded, t2, lb, ub)

	file_name = "result/anderson_tempo1_beta40.0_t$(t)_U$(U)_mu$(μ)_tep$(t2)_lb$(lb)_ub$(ub)_dt0.05_order$(order)_chi$(chi)_2.json"

	results = Dict("ns" => ns, "ts"=>τs, "gf"=>Gτ, "ws"=>ws, "Aw"=>Aw)

	println("save results to ", file_name)

	open(file_name, "w") do f
		write(f, JSON.json(results))
	end

	return results
end

function main2_all_U(t₀, t = t₀ + 20.; t2=100, lb=-2, ub=2, order=6, chi=1024)
	for U in [0.1, 0.5, 1.]
		main2(t, t₀; t2=t2, U=U, lb=lb, ub=ub, order=order, chi=chi)
	end
end