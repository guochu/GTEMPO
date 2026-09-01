# Fock-space representation of impurity Hamiltonians
# ---------------------------------------------------
# Jordan-Wigner matrices in the occupation number basis, and the mapping
# of a generic ImpurityHamiltonian to its Fock-space matrix. This is the
# input of the exact propagator construction in fockmatrix.jl
# (fock_propagator, fock_thermalstate).

const Id = [1 0; 0 1]
const Sz = [1 0; 0 -1]
const Sm = [0 1; 0 0] # Sm|1> = |0>
const Sp = [0 0; 1 0] # Sp|0> = |1>

const ket0 = [1; 0]
const ket1 = [0; 1]
const ket = (ket0, ket1)

function jw_operators(N::Int)
    annihilators = Vector{Matrix{Int}}(undef, N)
    creators = Vector{Matrix{Int}}(undef, N)

    # Sz ... Sz Sm/Sp Id ... Id
    for j in 1:N
        op_a = (j == 1) ? Sm : Sz
        op_adag = (j == 1) ? Sp : Sz

        for k in 2:N
            target_a = (k < j) ? Sz : (k == j ? Sm : Id)
            target_adag = (k < j) ? Sz : (k == j ? Sp : Id)

            op_a = kron(op_a, target_a)
            op_adag = kron(op_adag, target_adag)
        end

        annihilators[j] = op_a
        creators[j] = op_adag
    end

    return creators, annihilators
end

# fockstate(1,0,0,0) => band=1 occupied, band=1 is the most significant bit
function fockstate(ni::Int...)
    res = ket[ni[1] + 1]
    for i in 2:length(ni)
        res = kron(res, ket[ni[i] + 1])
    end
    return res
end

function fockmatrix(m::AndersonIM, bands::Int)
    adag, a = jw_operators(bands)
    if bands == 1
        return m.μ * adag[1]*a[1]
    elseif bands == 2
        return m.μ * (adag[1]*a[1] + adag[2]*a[2]) + m.U * adag[1]*a[1] * adag[2]*a[2]
    else
        error("Invalid bands of $bands")
    end
end

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
