
const Id = [1 0; 0 1]
const Sz = [1 0; 0 -1]
const Sm = [0 1; 0 0] # Sm|1> = |0>
const Sp = [0 0; 1 0] # Sm|0> = |1>

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

# fockstate(1,0,0,0) => band=1
# fockstate(0,1,0,0) => band=2
# fockstate(0,0,1,0) => band=3
# fockstate(0,0,0,1) => band=4
function fockstate(ni::Int...)
    res = ket[ni[1] + 1]
    for i in 2:length(ni)
        res = kron(res,ket[ni[i] + 1])
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
"""
	initfockstate2(lattice, fockstate::FockMatrix; normalize) -> GrassmannMPS

Convert the Fock-space matrix `fockstate` (an operator in the occupation
number basis, e.g. a density matrix) directly into a SparseGMPS (see
`_tosparsegmps`), which is multiplied into the vacuum state with the
sparse `mult!`. The bra variables sit on the forward branch and the ket
variables on the backward branch at time slice 1.
"""
function initfockstate2(lattice::RealGrassmannLattice, fockstate::FockMatrix; normalize::Bool=true)
    (fockstate.bands == lattice.bands) || throw(DimensionMismatch("FockMatrix bands $(fockstate.bands) do not match lattice bands $(lattice.bands)"))
    M = lattice.bands
    bpos = [index(lattice, 1, conj=true, branch=:+, band=i) for i in 1:M]
    kpos = [index(lattice, 1, conj=false, branch=:-, band=i) for i in 1:M]
    sparse = _tosparsegmps(lattice, fockstate, bpos, kpos)

    state = vacuumstate(lattice)
    if !(scalartype(sparse) <: Real) && (scalartype(state) <: Real)
        state = complex(state)
    end
    mult!(state, sparse, trunc=NoTruncation())

    alg = Orthogonalize(SVD(), normalize=normalize)
    open("/dev/null", "w") do devnull # slience the warning
        redirect_stderr(devnull) do
            canonicalize!(state, alg=alg)
        end
    end
    return state
end

"""
	initthermalstate2(lattice, model, β) -> GrassmannMPS

The impurity thermal equilibrium state `exp(-βĤ)/tr(exp(-βĤ))` (for
`β == Inf` the ground state projector) as a GrassmannMPS, built via
`fock_thermalstate` and `initfockstate2`.
"""
function initthermalstate2(lattice::RealGrassmannLattice, model, β::Real)
    initfockstate2(lattice, fock_thermalstate(model, β, lattice.bands); normalize=true)
end

