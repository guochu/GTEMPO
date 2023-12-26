function _g₁(β::Real, μ::Real, ε::Float64)
	if β == Inf
		if ε > μ
			return 1.0
		elseif ε < μ
			return 0.
		else
			return 0.5
		end
		# ifelse(ε >= μ, 1.0, 0.0)
	else
		1.0/(1.0+exp(-β * (ε-μ) ))
	end
end 
function _g₂(β::Real, μ::Real, ε::Float64)
	if β == Inf
		# ifelse(ε >= μ, 0.0, 1.0)
		if ε > μ
			return 0.0
		elseif ε < μ
			return 1.0
		else
			return 0.5
		end
	else
		1.0/(1.0+exp(β * (ε-μ) ))
	end
end 

const tol = 1.0e-6