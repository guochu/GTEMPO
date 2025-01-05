using ImpurityModelBase, JSON

J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(V::Real, t::Real, Lsys::Int; δt=0.1, dw=0.1)
	β = Inf
	D = 1.
	N = round(Int, t / δt)

	ts = [i*δt for i in 0:N]

	leftmu = V / 2
	rightmu = -V / 2


	leftbath_c = fermionicbath(spectrum_func(D), β=β, μ=leftmu)
	leftbath = discretebath(leftbath_c, δw=dw)
	rightbath_c = fermionicbath(spectrum_func(D), β=β, μ=rightmu)
	rightbath = discretebath(rightbath_c, δw=dw)

	J = 1
	hsys = zeros(Float64, Lsys, Lsys)
	for i in 1:Lsys-1
		hsys[i, i+1] = J
		hsys[i+1, i] = J
	end
	ρsys = zeros(ComplexF64, Lsys, Lsys)

	model = BoundaryDriving(hsys, leftbath, rightbath)
	h = cmatrix(model)
	ρ₀ = separablestate(model, ρsys)


	# leftcurrentop = leftparticlecurrent_cmatrix(model)
	# rightcurrentop = rightparticlecurrent_cmatrix(model)

	# cache = freefermions_cache(h)
	# n1 = ComplexF64[]
	# leftcurrent = ComplexF64[]
	# rightcurrent = ComplexF64[]

	# push!(n1, ρ₀[1,1])
	# push!(leftcurrent, sum(leftcurrentop .* ρ₀))
	# push!(rightcurrent, sum(rightcurrentop .* ρ₀))

	# for i in 2:length(ts)
	# 	ρ = freefermions_timeevo(ρ₀, h, ts[i], cache)
	# 	push!(n1, ρ[1,1])
	# 	push!(leftcurrent, sum(leftcurrentop .* ρ))
	# 	push!(rightcurrent, sum(rightcurrentop .* ρ))
	# end
	# return ts, n1, leftcurrent, rightcurrent

	gt, lt = freefermions_greater_lesser(h, ρ₀, ts, i=Lsys)
	return ts, gt
end