
function sysdynamics_util2(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model; idx::Int=1, branch::Symbol=:+, trunc::TruncationScheme=DefaultKTruncation)
    H = fockmatrix(model, lattice.bands)

	if branch == :+
	    coeff = - im * lattice.δt
	elseif branch == :-
	    coeff = im * lattice.δt
	elseif branch == :τ
	    coeff = - lattice.δτ
	end
	ib, ik = branch == :- ? (idx, idx+1) : (idx+1, idx)

    # direct exp will introduce many non-physcial non-zero elements
    # rho = exp(-β * H)
    vals, vecs = eigen(H)
    # E0 = minimum(vals)
    # exp_vals = exp.(coeff .* (vals .- E0))
    # exp_vals = exp_vals ./ sum(exp_vals)

	exp_vals = exp.(coeff .* vals)
    fockstate = vecs * diagm(exp_vals) * vecs'

    @assert lattice.bands == convert(Int, log2(size(fockstate,1)))
    terms = fock2grassmann(fockstate)

    vac = vacuumstate(lattice)
    states = map(terms) do term
        v, bind, kind = term.value, term.bra_indices, term.ket_indices
        bpos = [index(lattice, ib, conj=true, branch=branch, band=i) for i in bind]
        kpos = [index(lattice, ik, conj=false, branch=branch, band=i) for i in kind]
        apply!(GTerm(bpos..., kpos..., coeff=v), deepcopy(vac))
    end
    state = states[1]
    for i in 2:length(states)
        state = state + states[i]
    end

    alg = Orthogonalize(SVD(), normalize=false)
	canonicalize!(state, alg=alg)

	return mult(state, gmps, trunc=trunc)
end
function sysdynamics_forward2(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model; trunc::TruncationScheme=DefaultKTruncation)
	for i in 1:lattice.Nt
    	gmps = sysdynamics_util2(gmps, lattice, model, idx=i, branch=:+, trunc=trunc)
	end
	return gmps
end
function sysdynamics_backward2(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model; trunc::TruncationScheme=DefaultKTruncation)
	for i in 1:lattice.Nt
		gmps = sysdynamics_util2(gmps, lattice, model, idx=i, branch=:-, trunc=trunc)
	end
	return gmps
end
function sysdynamics_imaginary2(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model; trunc::TruncationScheme=DefaultKTruncation)
	for i in 1:lattice.Nτ
		gmps = sysdynamics_util2(gmps, lattice, model, idx=i, branch=:τ, trunc=trunc)
	end
	return gmps
end



function sysdynamics2(lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
    gmps = vacuumstate(lattice)
	return sysdynamics_imaginary2(gmps, lattice, model; trunc=trunc)
end 
function sysdynamics2(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
    gmps = vacuumstate(lattice)
    if isnothing(branch)
		gmps = sysdynamics_forward2(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward2(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward2(gmps, lattice, model; trunc=trunc) : sysdynamics_backward2(gmps, lattice, model; trunc=trunc)
	end
end 
function sysdynamics2(lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
    gmps = vacuumstate(lattice)
    if isnothing(branch)
		gmps = sysdynamics_forward2(gmps, lattice, model; trunc=trunc)
		gmps = sysdynamics_backward2(gmps, lattice, model; trunc=trunc)
		return sysdynamics_imaginary2(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward2(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward2(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary2(gmps, lattice, model; trunc=trunc)
		end
	end
end 



function sysdynamics_fast(lattice::ImagGrassmannLattice{O}, model::AbstractImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation) where O
	lattice2 = similar(lattice, ordering=A1B1B̄1Ā1())
	x = _sysdynamics_fast(lattice2, model, trunc=trunc)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end

_sysdynamics_fast(lattice::ImagGrassmannLattice{<:A1B1B̄1Ā1}, model::AbstractImpurityHamiltonian; kwargs...) = _sysdynamics_fast_timelocal2(lattice, model; kwargs...)

function _sysdynamics_fast_timelocal2(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; 
										trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice1 = similar(lattice, N=1)
	gmps1 = sysdynamics2(lattice1, model; trunc=trunc, kwargs...)
	posa = get_left(lattice1, 2)
	posb = get_right(lattice1, 1)
	function check_gmps1(x, tol)
		# println( space_l(x[posa]), " ", space_r(x[posb]))
		@assert isoneunit(space_l(x[posa]) )
		@assert isoneunit(space_r(x[posb]))
		s = 1.
		for i in 1:posa-1
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		for i in posb+1:length(x)
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		return s
	end
	s0 = check_gmps1(gmps1, 1.0e-10)
	_s = scaling(gmps1)^(length(gmps1) / (posb-posa+1))
	for i in posa:posb
		gmps1[i] *= _s
	end
	gmps1[posa] *= s0

	gmps = vacuumstate(lattice)
	for i in 1:lattice.N
		posa_n = get_left(lattice, i+1)
		posb_n = get_right(lattice, i)		
		for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx] 
		end
	end
	# return canonicalize!(gmps, alg=Orthogonalize(trunc=trunc))
	return gmps
end



function sysdynamics_fast(lattice::RealGrassmannLattice{O}, model::AbstractImpurityHamiltonian; 
							trunc::TruncationScheme=DefaultKTruncation, kwargs...) where {O}
	if LayoutStyle(lattice) isa TimeLocalLayout
		lattice2 = similar(lattice, ordering=A1B1ā1b̄1Ā1B̄1a1b1())
	else
		lattice2 = similar(lattice, ordering=A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2())
	end
	x = _sysdynamics_fast(lattice2, model; trunc=trunc, kwargs...)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end

_sysdynamics_fast(lattice::RealGrassmannLattice{<:A1B1ā1b̄1Ā1B̄1a1b1}, model::AbstractImpurityHamiltonian; kwargs...) = _sysdynamics_fast_timelocal2(lattice, model; kwargs...)

function _sysdynamics_fast(lattice::RealGrassmannLattice{<:A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}, model::AbstractImpurityHamiltonian; 
								trunc::TruncationScheme=DefaultKTruncation, kwargs...)
	lattice1 = similar(lattice, N=1)
	gmps1 = sysdynamics2(lattice1, model; trunc=trunc, kwargs...)
	posa_f = get_left(lattice1, 2, branch=:+)
	posb_f = get_right(lattice1, 1, branch=:+)
	posa_b = get_left(lattice1, 1, branch=:-)
	posb_b = get_right(lattice1, 2, branch=:-)
	function check_gmps1(x, tol)
		# println( space_l(x[posa_f]), " ", space_r(x[posb_f]), " ", space_l(x[posa_b]), " ", space_r(x[posb_b]))
		@assert isoneunit(space_l(x[posa_f]))
		@assert isoneunit(space_r(x[posb_f]))
		@assert isoneunit(space_l(x[posa_b]))
		@assert isoneunit(space_r(x[posb_b]))
		s = 1.
		for i in Iterators.flatten((1:posa_f-1, posb_f+1:posa_b-1, posb_b+1:length(x)))
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		return s
	end
	s0 = check_gmps1(gmps1, 1.0e-10)
	_s = scaling(gmps1)^(length(gmps1) / (posb_f-posa_f+1 + posb_b - posa_b+1))
	for i in Iterators.flatten((posa_f:posb_f, posa_b:posb_b))
		gmps1[i] *= _s
	end
	gmps1[posa_f] *= s0

	gmps = vacuumstate(lattice)
	for i in 1:lattice.N
		posa_n = get_left(lattice, i+1, branch=:+)
		posb_n = get_right(lattice, i, branch=:+)		
		for (idx, idx_n) in zip(posa_f:posb_f, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx] 
		end
		posa_n = get_left(lattice, i, branch=:-)
		posb_n = get_right(lattice, i+1, branch=:-)		
		for (idx, idx_n) in zip(posa_b:posb_b, posa_n:posb_n)
			gmps[idx_n] = gmps1[idx] 
		end		
	end

	return gmps
end


struct A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1 <: MixedGrassmannOrdering end
function index(x::MixedGrassmannLattice1Order{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ + 1) || throw(BoundsError(1:x.kτ, i))
			else
				(1 <= i <= x.Nt + 1) || throw(BoundsError(1:x.kt, i))
			end
		end
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*bands+1-band, band) 
	else
		k = x.Nt + 1
        if branch == :-
            ifelse(conj, 4i*bands+2*bands+band, 4i*bands+band) + (x.Nτ)*2*bands
        elseif branch == :+
            ifelse(conj, 4i*bands+bands+band, 4i*bands+3*bands+band) + (x.Nτ)*2*bands
        else
			ifelse(conj, (x.Nτ+1-i)*2*bands + 2bands+1-band, (x.Nτ+1-i)*2*bands + band) + 2*x.bands
		end
	end
end

get_left_imag(lattice::MixedGrassmannLattice{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, j::Int) = index(lattice, j, conj=true, band=lattice.bands, branch=:τ)
get_right_imag(lattice::MixedGrassmannLattice{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, j::Int) = index(lattice, j, conj=false, band=lattice.bands, branch=:τ)

get_left_real(lattice::MixedGrassmannLattice{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, j::Int) = index(lattice, j, conj=true, band=1, branch=:-)
get_right_real(lattice::MixedGrassmannLattice{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, j::Int) = index(lattice, j, conj=true, band=lattice.bands, branch=:+)


function sysdynamics_fast(lattice::MixedGrassmannLattice{O}, model::AbstractImpurityHamiltonian; 
							trunc::TruncationScheme=DefaultKTruncation, kwargs...) where {O}
	if LayoutStyle(lattice) isa TimeLocalLayout
		lattice2 = similar(lattice, ordering=A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1())
	else
		errro("sysdynamics_fast only be implemented for TimeLocalLayout MixedGrassmannLattice")
	end
	x = _sysdynamics_fast(lattice2, model; trunc=trunc, kwargs...)
	return changeordering(O, lattice2, x, trunc=trunc)[2]
end
function _sysdynamics_fast(lattice::MixedGrassmannLattice{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, model::AbstractImpurityHamiltonian; 
								trunc::TruncationScheme=DefaultKTruncation, branch::Union{Nothing, Symbol}=nothing)
    function check_gmps1(x, posa, posb, tol)
		# println( space_l(x[posa]), " ", space_r(x[posb]))
		@assert isoneunit(space_l(x[posa]) )
		@assert isoneunit(space_r(x[posb]))
		s = 1.
		for i in 1:posa-1
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		for i in posb+1:length(x)
			@assert isoneunit(space_l(x[i])) && isoneunit(space_r(x[i]))
			xi = only(x[i][(Z2Irrep(0), Z2Irrep(0), Z2Irrep(0))])
			s *= xi
			@assert abs(abs(xi) - 1) < tol
			# @assert abs(TK.scalar(x[i]) - 1) < tol
		end
		return s
	end
    function get_branch(branch::Symbol)
        lattice1 = similar(lattice, Nt=1, Nτ=1)
        gmps1 = sysdynamics2(lattice1, model; trunc=trunc, branch=branch)
        if branch == :τ
            posa, posb = get_left_imag(lattice1, 2), get_right_imag(lattice1, 1)
            N = lattice.Nτ
        else
            posa, posb = get_left_real(lattice1, 1), get_right_real(lattice1, 2)
            N = lattice.Nt
        end
        s0 = check_gmps1(gmps1, posa, posb, 1.0e-10)
        _s = scaling(gmps1)^(length(gmps1) / (posb-posa+1))
        for i in posa:posb
            gmps1[i] *= _s
        end
        gmps1[posa] *= s0

        gmps = vacuumstate(lattice)
        for i in 1:N
            if branch == :τ
                posa_n, posb_n = get_left_imag(lattice, i+1), get_right_imag(lattice, i)
            else
                posa_n, posb_n = get_left_real(lattice, i), get_right_real(lattice, i+1)
            end
            for (idx, idx_n) in zip(posa:posb, posa_n:posb_n)
                gmps[idx_n] = gmps1[idx] 
            end
        end
        return gmps
    end
    if isnothing(branch)
        gmpsp = get_branch(:+)
        gmpsm = get_branch(:-)
        gmpsτ = get_branch(:τ)
        gmps = mult(gmpsp, gmpsm, trunc=trunc)
        return mult(gmps, gmpsτ, trunc=trunc)
    else
        return get_branch(branch)
    end
end

