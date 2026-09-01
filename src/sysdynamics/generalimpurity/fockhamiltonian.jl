# Fock-space matrix of an ImpurityHamiltonian
# ------------------------------------------------
# The terms of a generic ImpurityHamiltonian are mapped to their
# Jordan-Wigner matrix representation; this is the input of the exact
# propagator construction in fockmatrix.jl (fock_propagator,
# fock_thermalstate).

function jw_operators(term::AdagATerm, N::Int)
	adag, a = jw_operators(N)
	return term.coeff * adag[term.positions[1]] * a[term.positions[2]]
end
function jw_operators(term::QuarticTerm, N::Int)
	adag, a = jw_operators(N)
	return term.coeff * adag[term.positions[1]] * adag[term.positions[2]] * a[term.positions[3]] * a[term.positions[4]]
end
function fockmatrix(h::ImpurityHamiltonian, bands::Int)
	@assert bands == h.bands
	adag, a = jw_operators(h.bands)
	sum(h.data) do data
		jw_operators(data, h.bands)
	end
end
