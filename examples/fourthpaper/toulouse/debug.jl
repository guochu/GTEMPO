using JSON, QuadGK


function J(D::Real, ω::Real)
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

function Gw(ω::Float64, spectrum, ϵ_d, μ=0)
	δ = 1e-8
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
	1.0/(ω-ϵ_d-quadgk(ε -> f(ε)/(ω+μ-ε+im*δ), lb, ub)[1])
end

function G(t::Float64, spectrum, ϵ_d, μ=0)
    δ = 1e-8
    A = quadgk(ω -> (Gw(ω, spectrum, ϵ_d, μ)-1.0/(ω+im*δ))*exp(-im*ω*t), -20., 20.)[1]
    im*(A/(2π)-im)
end

function main_analytic(t::Float64; ϵ_d=0, δt=0.05)
	D = 2.
	N = round(Int, t / δt)
	ts = [i*δt for i in 0:N]
	g = [G(t, spectrum_func(D), -ϵ_d, 0.) for t in ts]	

	g = real(g)
	data_path = "result/analytic_mu$(ϵ_d)_dt$(δt)_N$(N).json"

	results = Dict("ts"=>ts, "gf" => g)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g
end