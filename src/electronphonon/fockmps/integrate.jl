
"""
	integrate(x::FockMPS)

Fully contract the state represented by the MPS and return the resulting overall scalar value (equivalent to the total coefficient obtained by summing over all sites).
Delegates to `sum` of the inner FiniteMPSAlgorithms `CanonicalMPS` payload (including the `scaling` factor).
"""
integrate(x::FockMPS) = sum(x.parent)
