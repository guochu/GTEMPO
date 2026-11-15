# 决定性实验:U1GTEMPO2 参照 K' 链(奇插入)喂 GM.Zvalue
using Pkg
Pkg.activate("/home/guochu/Documents/Missile/U1GTEMPO2")
using U1GTEMPO2
import GrassmannMPS as GM
using LinearAlgebra: Diagonal, tr
μ = 0.7; β = 1.0
q = exp(-β * μ)
p0 = [1.0, q, q, q * q]; p0 ./= sum(p0)
c = RealTime(Nt=4, dt=0.15)
bath = DiscreteBath([1.0], [0.5]; beta=β)
p = SIAMProblem(c, SIAM(energy=(μ, μ), U=0.0); bath=(bath, bath), order=2,
                initial_probabilities=Tuple(p0), alg=PartialIF())
s = solve(p)
cn = 1.0
o = U1GTEMPO2.fock_operators()
ins = U1GTEMPO2.correlation_insertions(:+, 1, :+, 0, :greater, :up)
println("插入 = ", [(i.branch, i.boundary, i.charge) for i in ins])
Kref = U1GTEMPO2.system_factor_reference(c, SIAM(energy=(μ, μ), U=0.0); order=2,
        initial_probabilities=Tuple(p0), coupling_norm=cn, insertions=ins)
Is = Tuple(f.state for f in p.influence)
z = GM.Zvalue([Kref, Is...])
Z = s.Z
ratio = exp(z.logabs) * z.phase / Z
println("GM.Zvalue(Kref, IFs)/Z = ", ratio, "   logabs = ", z.logabs)
println("参照 greater(1)/(-im)  = ", greater(s, 1; spin=:up) / (-1.0im))
println("若两者一致 → GM 引擎无罪,port pairize 为凶;若 GC 值错 → GM 引擎为凶")
