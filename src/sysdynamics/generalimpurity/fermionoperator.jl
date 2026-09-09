# abstract type AbstractFTerm end

positions(x::AbstractTerm) = x.positions

# """
# 	struct TunnelingTerm <: AbstractFTerm

# Fermionic twobody term
# """
# struct TunnelingTerm{T <: Number} <: AbstractFTerm
# 	positions::Tuple{Int, Int}
# 	coeff::T
# end

# TunnelingTerm(pos::Tuple{Int, Int}; coeff::Number=1) = TunnelingTerm(pos, float(coeff))
# TunnelingTerm(i::Int, j::Int; kwargs...) = TunnelingTerm((i, j); kwargs...)

# tunneling(i::Int, j::Int; kwargs...) = TunnelingTerm(i, j; kwargs...)

# TK.scalartype(::Type{TunnelingTerm{T}}) where {T} = T

# function Base.adjoint(x::TunnelingTerm)
# 	i, j = positions(x)
# 	return TunnelingTerm((j, i), coeff=conj(x.coeff))
# end

# Base.copy(x::TunnelingTerm) = TunnelingTerm(positions(x), copy(x.coeff))

# Base.:*(s::TunnelingTerm, m::Number) = TunnelingTerm(positions(s), coeff=s.coeff * m)
# Base.:*(m::Number, s::AbstractFTerm) = s * m
# Base.:/(s::AbstractFTerm, m::Number) = s * (1 / m)
# Base.:+(s::AbstractFTerm) = s
# Base.:-(s::AbstractFTerm) = (-1) * s


# """
# 	struct FourBodyTerm <: AbstractFTerm

# Fermionic fourbody term, ĉ₁†ĉ₂†ĉ₃ĉ₄
# """
# struct InteractionTerm{T<:Number} <: AbstractFTerm
# 	positions::NTuple{4, Int}
# 	coeff::T
# end

# InteractionTerm(pos::NTuple{4, Int}; coeff::Real=1) = InteractionTerm(pos, float(coeff))
# InteractionTerm(i::Int, j::Int, k::Int, l::Int; kwargs...) = InteractionTerm((i, j, k, l); kwargs...)

# interaction(i::Int, j::Int, k::Int, l::Int; kwargs...) = InteractionTerm(i, j, k, l; kwargs...)

# TK.scalartype(::Type{InteractionTerm{T}}) where {T} = T

# function Base.adjoint(x::InteractionTerm)
# 	i, j, k, l = positions(x)
# 	return interaction((l,k,j,i), coeff=conj(x.coeff))
# end

# Base.:*(s::InteractionTerm, m::Number) = InteractionTerm(positions(s), coeff=s.coeff * m)


"""
	struct ImpurityHamiltonian <: ConstImpurityHamiltonian

Generic constant impurity Hamiltonian built from a list of `NormalTerm`s
(`tunneling` / `interaction`), with `bands` impurity bands. Terms can be
appended with `push!` or set through the constructor.
"""
struct ImpurityHamiltonian <: ConstImpurityHamiltonian
	data::Vector{NormalTerm}
	bands::Int

function ImpurityHamiltonian(data::Vector{<:NormalTerm}, bands::Int)
	@boundscheck begin
		for f in data
			for j in positions(f)
				(1 <= j <= bands) || throw(BoundsError(1:bands, j))
			end
		end
	end
	new(convert(Vector{NormalTerm}, data), bands)
end

end
ImpurityHamiltonian(data::Vector{<:NormalTerm}; bands::Int=1) = ImpurityHamiltonian(data, bands)
ImpurityHamiltonian(; bands::Int = 1) = ImpurityHamiltonian(Vector{NormalTerm}(), bands)
function Base.push!(x::ImpurityHamiltonian, f::NormalTerm)
	@boundscheck begin
		for j in f.positions
			(1 <= j <= x.bands) || throw(BoundsError(1:x.bands, j))
		end
	end
	# for t in x.data
	# 	if positions(t) == positions(f)
	# 		throw(ArgumentError("duplicate positions not allowed"))
	# 	end
	# end
	push!(x.data, f)
end 

function TK.scalartype(h::ImpurityHamiltonian)
	T = Float64
	for f in h.data
		T = promote_type(T, eltype(f))
	end
	return T
end

"""
	struct QuenchedImpurityHamiltonian <: ConstImpurityHamiltonian

Impurity Hamiltonian with a quench protocol: the impurity evolves with the
pre-quench Hamiltonian `hτ` on the imaginary-time branch (`:τ`, and for the
thermal initial state) and with the post-quench Hamiltonian `ht` on the
real-time branches (`:+`/`:-`). The terms can only be set through the
constructor, there is no `push!`.
"""
struct QuenchedImpurityHamiltonian <: ConstImpurityHamiltonian
	hτ::Vector{NormalTerm}
	ht::Vector{NormalTerm}
	bands::Int

	function QuenchedImpurityHamiltonian(hτ::Vector{<:NormalTerm}, ht::Vector{<:NormalTerm}, bands::Int)
		@boundscheck begin
			for f in hτ
				for j in positions(f)
					(1 <= j <= bands) || throw(BoundsError(1:bands, j))
				end
			end
			for f in ht
				for j in positions(f)
					(1 <= j <= bands) || throw(BoundsError(1:bands, j))
				end
			end
		end
		new(convert(Vector{NormalTerm}, hτ), convert(Vector{NormalTerm}, ht), bands)
	end

end
QuenchedImpurityHamiltonian(hτ::Vector{<:NormalTerm}, ht::Vector{<:NormalTerm}; bands::Int=1) =
	QuenchedImpurityHamiltonian(hτ, ht, bands)

function TK.scalartype(h::QuenchedImpurityHamiltonian)
	T = Float64
	for f in h.hτ
		T = promote_type(T, eltype(f))
	end
	for f in h.ht
		T = promote_type(T, eltype(f))
	end
	return T
end

"""
	struct TdImpurityOp

Time-dependent collection of `NormalTerm`s: `data` multiplied by the scalar
`f(t)` gives the `ImpurityHamiltonian` contribution of this term at time `t`.
"""
struct TdImpurityOp
	data::Vector{NormalTerm}
	f::Function
	bands::Int

	function TdImpurityOp(data::Vector{<:NormalTerm}, f::Function, bands::Int)
		@boundscheck begin
			for term in data
				for j in positions(term)
					(1 <= j <= bands) || throw(BoundsError(1:bands, j))
				end
			end
		end
		new(convert(Vector{NormalTerm}, data), f, bands)
	end

end
TdImpurityOp(data::Vector{<:NormalTerm}, f::Function; bands::Int=1) = TdImpurityOp(data, f, bands)

function (x::TdImpurityOp)(t::Real)
	c = x.f(t)
	data = [item * c for item in x.data]
	return ImpurityHamiltonian(data, x.bands)
end

# hτ	is used for the imaginary-time branch of the contour
# ht	is used for the constant part of the real-time branch of the contour
# htt	is used for the time-dependent part of the real-time branch of the contour
"""
	struct TdImpurityHamiltonian <: AbstractTdImpurityHamiltonian

Time-dependent impurity Hamiltonian. `hτ` is used for the imaginary-time
branch of the contour (and for the thermal initial state), while the
real-time branches (`:+`/`:-`) evolve with the time-dependent Hamiltonian
`ht + Σ htt(t)`: calling the model with a time `t` returns the constant
`ImpurityHamiltonian` of the real-time branches at `t`.
"""
struct TdImpurityHamiltonian <: AbstractTdImpurityHamiltonian
	hτ::Vector{NormalTerm}
	ht::Vector{NormalTerm}
	htt::Vector{TdImpurityOp}
	bands::Int

	function TdImpurityHamiltonian(hτ::Vector{<:NormalTerm}, ht::Vector{<:NormalTerm},
									htt::Vector{<:TdImpurityOp}, bands::Int)
		@boundscheck begin
			for f in (hτ..., ht...)
				for j in positions(f)
					(1 <= j <= bands) || throw(BoundsError(1:bands, j))
				end
			end
			for op in htt
				(op.bands == bands) || throw(DimensionMismatch("TdImpurityOp bands $(op.bands) do not match model bands $bands"))
			end
		end
		new(convert(Vector{NormalTerm}, hτ), convert(Vector{NormalTerm}, ht),
			convert(Vector{TdImpurityOp}, htt), bands)
	end

end
TdImpurityHamiltonian(hτ::Vector{<:NormalTerm}, ht::Vector{<:NormalTerm},
						htt::Vector{<:TdImpurityOp}; bands::Int=1) =
	TdImpurityHamiltonian(hτ, ht, htt, bands)

function (x::TdImpurityHamiltonian)(t::Real)
	ht = copy(x.ht)
	for op in x.htt
		append!(ht, op(t).data)
	end
	return ImpurityHamiltonian(ht, x.bands)
end