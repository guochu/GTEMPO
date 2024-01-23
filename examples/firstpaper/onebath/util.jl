push!(LOAD_PATH, "../../../src")

using GTEMPO


function J(D::Real, ω::Real)
	t′ = 0.3162 
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * t′^2 
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)
