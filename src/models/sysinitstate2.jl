
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
function fockmatrix(m::KanamoriIM, bands::Int)
    norb, U, J, μ = m.norb, m.U, m.J, m.μ
	@assert bands == 2*norb
    N = 2 * norb
    adag, a = jw_operators(N)

    mat = μ * sum(adag[i]*a[i] for i in 1:N)
    for x in 1:m.norb
        xu, xd = 2*x-1, 2*x
        mat += U * adag[xu]*adag[xd]*a[xd]*a[xu]
        for y in 1:m.norb
            yu, yd = 2*y-1, 2*y
            if x != y
                mat += (U - 2*J) * adag[xu]*adag[yd]*a[yd]*a[xu]
                mat -= J * (adag[xu]*adag[xd]*a[yu]*a[yd] + adag[xu]*adag[yd]*a[yu]*a[xd])
            end
            if x > y
                mat += (U - 3*J) * adag[xu]*adag[yu]*a[yu]*a[xu]
                mat += (U - 3*J) * adag[xd]*adag[yd]*a[yd]*a[xd]
            end
        end
    end
    return mat
end

function fock2grassmann(H::AbstractMatrix; δ::Float64 = 1e-10)
    M::Int = convert(Int, log2(size(H,1)))
    coherent_terms = []
    for row in 1:2^M
        for col in 1:2^M
            val = H[row, col]
            abs(val) < δ && continue
            m_int = row - 1
            n_int = col - 1
            
            m_bits = reverse(digits(m_int, base=2, pad=M))
            n_bits = reverse(digits(n_int, base=2, pad=M))
            
            bra_indices = [i for i in 1:M if m_bits[i] == 1] # order of 1 -> M
            ket_indices = [i for i in reverse(1:M) if n_bits[i] == 1]
            if length(bra_indices) != length(ket_indices)
                @warn "Ignore non-physical element: $bra_indices $ket_indices => $val"
                continue
            end
            push!(coherent_terms, (
                value = val, 
                bra_indices = bra_indices, 
                ket_indices = ket_indices
            ))
        end
    end
    return coherent_terms
end


function initfockstate(lattice::RealGrassmannLattice, fockstate::AbstractMatrix; normalize::Bool=true)
    @assert lattice.bands == convert(Int, log2(size(fockstate,1)))
    terms = fock2grassmann(fockstate)

    vac = vacuumstate(lattice)
    states = map(terms) do term
        v, bind, kind = term.value, term.bra_indices, term.ket_indices
        bpos = [index(lattice, 1, conj=true, branch=:+, band=i) for i in bind]
        kpos = [index(lattice, 1, conj=false, branch=:-, band=i) for i in kind]
        apply!(GTerm(bpos..., kpos..., coeff=v), deepcopy(vac))
    end
    state = states[1]
    for i in 2:length(states)
        state = state + states[i]
    end

    alg = Orthogonalize(SVD(), normalize=normalize)
    open("/dev/null", "w") do devnull # slience the warning
        redirect_stderr(devnull) do
            canonicalize!(state, alg=alg)
        end
    end
    return state
end

function initthermalstate(lattice::RealGrassmannLattice, model, β::Real)
    H = fockmatrix(model, lattice.bands)

    # direct exp will introduce many non-physcial non-zero elements
    # rho = exp(-β * H)
    vals, vecs = eigen(H)
	if β == Inf
		rho = vecs[:,1] * vecs[:,1]'
	else
		E0 = minimum(vals)
		exp_vals = exp.(-β .* (vals .- E0))
		exp_vals = exp_vals ./ sum(exp_vals)
		rho = vecs * diagm(exp_vals) * vecs'
	end

    initfockstate(lattice, rho; normalize=true)
end

