push!(LOAD_PATH, "../../../src")
using GTEMPO



J(ε) = ε^3*exp(-ε)

function main_imag_analytic(ϵ_d; β=1, N=20)
	spec = SpectrumFunction(J, lb = 0, ub =Inf)
	g = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, Nτ=N)
	return g
end