# insert an n̂ operator at time step i into the bare impurity dynamics, which is used to measure density-density correlations


function insert_n!(lattice::AbstractGrassmannLattice, gmps::GrassmannMPS, j::Int; band::Int=1, branch=:τ)
	@assert j > 0
	pos = index(lattice, j, branch=branch, conj=false, band=band)
 	t = copy(gmps[pos])
 	for (fl, fr) in fusiontrees(t)
 		if iseven(fl.uncoupled[2].n)
 			fill!(t[fl, fr], 0)
 		end
 	end
 	gmps[pos] = t
 	unset_svectors!(gmps)
 	return gmps
end 

insert_n(lattice::AbstractGrassmannLattice, gmps::GrassmannMPS, j::Int; kwargs...) = insert_n!(lattice, copy(gmps), j; kwargs...)