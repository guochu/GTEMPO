"""
	mult(x::FockMPS, y::FockMPS; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)

Compute the (compressed) product of two FockMPSs and return a new `FockMPS` (inputs are not modified).
The SVD compression route delegates to FiniteMPSAlgorithms' `hadamard`: the exact
(Hadamard) product `x ⊙ y` is formed and compressed by a single SVD sweep under `trunc`.

# Returns
The product `FockMPS`.
"""
mult(x::FockMPS, y::FockMPS; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0) =
	FockMPS(FMA.hadamard(x.parent, y.parent, FMA.SVDCompression(trunc=_fmatrunc(trunc), verbosity=verbosity)))

"""
	mult!(x::FockMPS, y::FockMPS; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)

Compute the (compressed) product of two FockMPSs in place on `x`, and return `x`.
"""
mult!(x::FockMPS, y::FockMPS; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0) =
	copy!(x, mult(x, y; trunc, verbosity))
