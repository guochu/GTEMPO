# termdynamics(lattice::AbstractGrassmannLattice, h::AbstractFTerm; kwargs...) = termdynamics!(vacuumstate(lattice), lattice, h; kwargs...)

# """
# 	termdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::AbstractFTerm; trunc)
# """
# function termdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, h::AbstractFTerm; kwargs...)
# 	for timestep in 1:lattice.Nτ
# 		termdynamics!(gmps, lattice, h, timestep; kwargs...)
# 	end
# 	return gmps
# end
# function termdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, h::AbstractFTerm, timestep::Int; trunc::TruncationScheme=DefaultKTruncation)
# 	all(x->1<=x<=lattice.bands, positions(h)) || throw(BoundsError(1:lattice.bands, positions(h)))
# 	(1 <= timestep <= lattice.Nτ) || throw(BoundsError(1:lattice.Nτ, timestep))
# 	return termdynamics_util!(gmps, lattice, h, timestep, -lattice.δτ, :τ, trunc)
# end

# function termdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, h::AbstractFTerm; 
# 						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
# 	if isnothing(branch)
# 		for timestep in 1:lattice.Nt
# 			termdynamics!(gmps, lattice, h, timestep; branch=:+, trunc=trunc)
# 			termdynamics!(gmps, lattice, h, timestep; branch=:-, trunc=trunc)
# 		end
# 	else
# 		for timestep in 1:lattice.Nt
# 			termdynamics!(gmps, lattice, h, timestep; branch=branch, trunc=trunc)
# 		end
# 	end
# 	return gmps
# end
# function termdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, h::AbstractFTerm, timestep::Int; 
# 						branch::Symbol, trunc::TruncationScheme=DefaultKTruncation)
# 	all(x->1<=x<=lattice.bands, positions(h)) || throw(BoundsError(1:lattice.bands, positions(h)))
# 	(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
# 	(1 <= timestep <= lattice.Nt) || throw(BoundsError(1:lattice.Nt, timestep))
# 	if branch == :+
# 		return termdynamics_util!(gmps, lattice, h, timestep, -im*lattice.δt, :+, trunc)
# 	else
# 		return termdynamics_util!(gmps, lattice, h, timestep, im*lattice.δt, :-, trunc)
# 	end
# end

# function termdynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, h::AbstractFTerm; 
# 						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
# 	if isnothing(branch)
# 		for timestep in 1:lattice.Nτ
# 			termdynamics!(gmps, lattice, h, timestep; branch=:τ, trunc=trunc)
# 		end
# 		for timestep in 1:lattice.Nt
# 			termdynamics!(gmps, lattice, h, timestep; branch=:+, trunc=trunc)
# 			termdynamics!(gmps, lattice, h, timestep; branch=:-, trunc=trunc)
# 		end
# 	else
# 		if branch == :τ
# 			for timestep in 1:lattice.Nτ
# 				termdynamics!(gmps, lattice, h, timestep; branch=branch, trunc=trunc)
# 			end
# 		else
# 			for timestep in 1:lattice.Nt
# 				termdynamics!(gmps, lattice, h, timestep; branch=branch, trunc=trunc)
# 			end
# 		end
# 	end
# 	return gmps
# end
# function termdynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, h::AbstractFTerm, timestep::Int; 
# 						branch::Symbol, trunc::TruncationScheme=DefaultKTruncation)
# 	all(x->1<=x<=lattice.bands, positions(h)) || throw(BoundsError(1:lattice.bands, positions(h)))
# 	(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
# 	if branch == :+
# 		(1 <= timestep <= lattice.Nt) || throw(BoundsError(1:lattice.Nt, timestep))
# 		return termdynamics_util!(gmps, lattice, h, timestep, -im*lattice.δt, :+, trunc)
# 	elseif branch == :-
# 		(1 <= timestep <= lattice.Nt) || throw(BoundsError(1:lattice.Nt, timestep))
# 		return termdynamics_util!(gmps, lattice, h, timestep, im*lattice.δt, :-, trunc)
# 	else
# 		(1 <= timestep <= lattice.Nτ) || throw(BoundsError(1:lattice.Nτ, timestep))
# 		return termdynamics_util!(gmps, lattice, h, timestep, -lattice.δτ, :τ, trunc)
# 	end
# end


# function termdynamics_util!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::AbstractFTerm, j::Int, dt::Number, branch::Symbol, trunc)
# 	alg = Orthogonalize(trunc = trunc)
# 	m = get_Gterm(lattice, h, j, dt, branch)
# 	apply!(m, gmps)
# 	canonicalize!(gmps, alg=alg)
# 	return gmps
# end