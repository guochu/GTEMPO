# port builder 忠实复刻:含 cut-0 折叠 + prefix 机制
import GTEMPO
using GTEMPO
using GTEMPO: sectortype
using ImpurityModelBase
import GrassmannMPS as GM
using LinearAlgebra: Diagonal, I, tr, norm
μ = 0.7; β = 1.0
bath = discretefermionicbath([1.0], [0.5]; β=β)
lat = GrassmannLattice(N=4, δt=0.15, bands=2, contour=:real, symmetry=U1Symmetry())
q = exp(-β * μ)
ρ0 = Diagonal(ComplexF64[1.0, q, q, q * q]); ρ0 /= tr(ρ0)
spec = GTEMPO.sectortype(U1Symmetry())
qc = GTEMPO.var_charge(spec, 1, 2)
d = GTEMPO._fock_ops(2)[1]; ddag = d'
weights = GTEMPO.local_weights(spec, 2, 2)
fqs = GTEMPO.fock_charges(spec, 2)
prs = [(s, t) for t in 1:4 for s in 1:4]
cpc = GTEMPO.cell_charges(spec, 2)
model = AndersonIM(U=0.0, ϵ_d=μ)

function faithful(insertions, tag)
    cuts = Dict{Int,Vector{OperatorInsertion}}()
    for ins in insertions; push!(get!(cuts, GTEMPO.insertion_cut_index(lat, ins), OperatorInsertion[]), ins); end
    cutleft(k) = foldr((ins, A) -> ins.matrix * A, filter(i -> i.side === :left, get(cuts, k, OperatorInsertion[])); init=Matrix{ComplexF64}(I, 4, 4))
    cutright(k) = foldl((A, ins) -> A * ins.matrix, filter(i -> i.side === :right, get(cuts, k, OperatorInsertion[])); init=Matrix{ComplexF64}(I, 4, 4))
    rho = cutleft(0) * Matrix{ComplexF64}(ρ0) * cutright(0)
    zq = GTEMPO.zero_charge(spec)
    prefix = zq
    cells = GTEMPO.chronological_cells(lat)
    L = length(cells)
    data = nothing
    for (k, (b, j)) in enumerate(cells)
        dt = lat.δt
        t = b === :+ ? (j - 1) * dt : (lat.Nt - j) * dt
        S = GTEMPO._model_propagator(model, b, dt / 2, t)
        mats = [S * M * S for M in weights]
        Ol = cutleft(k); Or = cutright(k)
        mats = [Ol * M * Or for M in mats]
        oddlater = 0
        for (cut, items) in cuts, ins in items
            cut >= k || continue
            oddlater += GTEMPO.parity_label(ins.charge === nothing ? zq : ins.charge) ? 1 : 0
        end
        for m in 0:15; mats[m+1] .*= (-1)^(oddlater * count_ones(m)); end
        leftlabels = k == 1 ? [zq] : [GTEMPO.qsub(GTEMPO.qsub(fqs[s], fqs[t]), prefix) for (s, t) in prs]
        prefix = zq
        for (cut, items) in cuts, ins in items
            cut == k || continue
            ins.charge === nothing || (prefix = GTEMPO.qadd(prefix, ins.charge))
        end
        rightlabels = k == L ? [zq] : [GTEMPO.qsub(GTEMPO.qsub(fqs[s], fqs[t]), prefix) for (s, t) in prs]
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
        println("  $tag cell k=$k: nnz(>0.05)=", count(abs.(core) .> 0.05), "  prefix=", prefix)
        ts = GTEMPO.split_cell(spec, core, leftlabels, rightlabels, cpc)
        data = data === nothing ? ts : vcat(data, ts)
    end
    out = GTEMPO.pairize(spec, data)
    GM.reorder!(out, collect(1:length(out)); trunc=GM.truncrelerr(1e-14))
    println("  $tag: logscale = ", round(GM.logscale(out), digits=2))
    return out
end

# A) U1GTEMPO2 模式:ddag@cut0(+q,折 ρ₀),d@(−,1)(−q)
fA = faithful([OperatorInsertion(:+, 0, ddag; charge=qc), OperatorInsertion(:-, 1, d; charge=GTEMPO.qneg(qc))], "U1模式")
# B) port 模式:ddag@cut0(−q,折 ρ₀),d@(+,1)(+q)
fB = faithful([OperatorInsertion(:+, 0, ddag; charge=GTEMPO.qneg(qc)), OperatorInsertion(:+, 1, d; charge=qc)], "port模式")
# C) 单插入:d@(−,1)(−q),无 cut0
fC = faithful([OperatorInsertion(:-, 1, d; charge=GTEMPO.qneg(qc))], "单d(−,1)")
println("三种模式链均构建完成")

# 积分忠实链:Z(fA·IFs)/Z(K·IFs) 应 = greater/(-im) = 0.66081-0.06945i
import GTEMPO: _inserted_expectation
Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=1.0), 2)
K = sysdynamics_weyl2order(lat, AndersonIM(U=0.0, ϵ_d=μ); coupling_norm=1.0, initial_state=ρ0)
zs = integrate_scaled(K, Is...); Z = exp(zs.logabs) * zs.phase
for (f, tag) in ((fA, "U1模式"), (fB, "port模式"), (fC, "单d(−,1)"))
    z = integrate_scaled(f, Is...)
    println(tag, ": Z'/Z = ", exp(z.logabs) * z.phase / Z)
end
println("参照: greater(1)/(-im) = ", (-0.06945314205699218 - 0.6608068990426897im) / (-1.0im))

# 迹:Integrate faithful U1-pattern chain + IFs vs 基准
Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=1.0), 2)
K = sysdynamics_weyl2order(lat, AndersonIM(U=0.0, ϵ_d=μ); coupling_norm=1.0, initial_state=ρ0)
zs = integrate_scaled(K, Is...); Z = exp(zs.logabs) * zs.phase
zf = integrate_scaled(fA, Is...)
println("U1 faithful: Z'/Z = ", exp(zf.logabs) * zf.phase / Z)
println("参照值 = 0.6608068990426886 - 0.06945314205699178im")
