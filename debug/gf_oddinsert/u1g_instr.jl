# U1GTEMPO2 system_factor_reference 逐胞仪器化(复刻其源码循环)
using Pkg
Pkg.activate("/home/guochu/Documents/Missile/U1GTEMPO2")
using U1GTEMPO2, ImpurityModelBase
import TensorKit as TK
using LinearAlgebra: Diagonal, I, tr

μ = 0.7; β = 1.0
q = exp(-β * μ)
c = RealTime(Nt=4, dt=0.15)
h = SIAM(energy=(μ, μ), U=0.0)
p0 = [1.0, q, q, q * q]; p0 ./= sum(p0)
order = 2; cn = 1.0
o = U1GTEMPO2.fock_operators(); d_up = o.up; d_up_dag = d_up'
insertions = [U1GTEMPO2.OperatorInsertion(:+, 0, d_up_dag; charge=(1, 0)),
              U1GTEMPO2.OperatorInsertion(:-, 1, d_up; charge=(-1, 0))]
cuts = Dict{Int,Vector{U1GTEMPO2.OperatorInsertion}}()
for ins in insertions; push!(get!(cuts, U1GTEMPO2.insertion_cut(c, ins), U1GTEMPO2.OperatorInsertion[]), ins); end
reduce(U1GTEMPO2.qadd, (ins.charge for ins in insertions); init=(0, 0)) == (0, 0) || error("charge")
rho = Matrix(Diagonal(ComplexF64.(p0)))
cutop(k) = foldr((ins, A) -> ins.matrix * A, get(cuts, k, U1GTEMPO2.OperatorInsertion[]); init=Matrix{ComplexF64}(I, 4, 4))
rho = cutop(0) * rho
weights = U1GTEMPO2.local_weight(order)
fqs = [(0, 0), (1, 0), (0, 1), (1, 1)]
prs = [(s, t) for t in 1:4 for s in 1:4]
cells = U1GTEMPO2.chronological_cells(c)
L = length(cells)
function build_chain()
    global data = nothing
    prefix = (0, 0)
    data = nothing
    for (k, (b, j)) in enumerate(cells)
        lam = U1GTEMPO2.step(c, b); H = U1GTEMPO2.hamiltonian(h)
        S = exp(-lam * H / (order == 2 ? 2 : 1))
        mats = order == 2 ? [S * M * S for M in weights] : b == :- ? [M * S for M in weights] : [S * M for M in weights]
        O = cutop(k); mats = [O * M for M in mats]
        oddlater = sum((isodd(sum(ins.charge)) for (cut, items) in cuts if cut >= k for ins in items); init=0)
        for m in 0:15; mats[m+1] .*= (-1)^(oddlater * count_ones(m)); end
        leftlabels = k == 1 ? [(0, 0)] : [U1GTEMPO2.qsub(U1GTEMPO2.qsub(fqs[s], fqs[t]), prefix) for (s, t) in prs]
        prefix = reduce(U1GTEMPO2.qadd, (ins.charge for ins in get(cuts, k, U1GTEMPO2.OperatorInsertion[])); init=prefix)
        rightlabels = k == L ? [(0, 0)] : [U1GTEMPO2.qsub(U1GTEMPO2.qsub(fqs[s], fqs[t]), prefix) for (s, t) in prs]
        core = zeros(ComplexF64, length(leftlabels), 16, length(rightlabels))
        for m in 0:15
            M = mats[m+1]
            if k == 1
                for (r, (s, t)) in enumerate(prs); core[1, m+1, r] = (M * rho)[s, t]; end
            elseif k == L
                for (l, (s, t)) in enumerate(prs); core[l, m+1, 1] = M[t, s]; end
            else
                for (l, (s, t)) in enumerate(prs), (r, (v, w)) in enumerate(prs)
                    t == w && (core[l, m+1, r] = M[v, s])
                end
            end
        end
        nz = count(abs.(core) .> 0.05)
        println("U1 cell k=$k ($b,$j): prefix=$prefix  nnz(>0.05)=$nz  nleft=$(length(leftlabels)) nright=$(length(rightlabels))")
        ts = U1GTEMPO2.split_cell(core, leftlabels, rightlabels, [(1, 0), (-1, 0), (0, 1), (0, -1)])
        data = data === nothing ? ts : vcat(data, ts)
    end
    println("U1 总变量张量数 = ", length(data))
end
build_chain()
