push!(LOAD_PATH, "../../../src")

using ImpurityModelBase
using JSON, Interpolations


function main_analytic(ϵ_d; ω₀=1, α₀=0.5)
	lb = -20
	ub = 20
	dw = 1.0e-2

	spec = semicircular(t=1)

	freqs = collect(frequencies(lb=lb, ub=ub, dw=dw))
	Gw = [holstein_Gw(spec, ϵ, g=-sqrt(α₀), ω=ω₀, ϵ_d=ϵ_d, order=20) for ϵ in freqs]
	println("Gw[1]=", Gw[1], ", Gw[end]=", Gw[end])

	return freqs, Gw
end

parse_complex_number(a) = a["re"] + im * a["im"]

function main(ϵ_d; β=Inf, t=1, N=100, ω₀=1, α₀=0.5, chi = 80)
	δt = t / N
	file_name = "result/holstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_mu$(ϵ_d)_chi$(chi).json"
	data = JSON.parsefile(file_name)

	ts = [data["ts"]...]
	gt = [-im*parse_complex_number(item) for item in data["gt"]]

	gt′ = linear_predict(gt, δt, tol=1.0e-9)

	println(gt′[end-5:end])

	t_all = δt * (length(gt′)-1)
	println("tall is ", t_all)

	interp=linear_interpolation(0:δt:t_all, gt′)
	δt′ = 0.001
	ts′ = 0:δt′:t_all
	gf′ = interp.(ts′)

	# return ts, gt, ts′, gf′

	lb = -20
	ub = 20
	dw = 1.0e-2
	Gw = Gt_to_Gw(gf′, δt′, lb=lb, ub=ub, dw=dw)
	println("Gw[1]=", Gw[1], ", Gw[end]=", Gw[end])

	freqs = collect(frequencies(lb=lb, ub=ub, dw=dw))

	return freqs, Gw
end



