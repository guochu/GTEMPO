# ============================================================================
# 性能对比基准（Toulouse 模型 / 非相互作用 SIAM）：PartialIF vs XTRGIF vs ExactTTIIF vs TDVPIF
#
# 模型：Toulouse 模型（ϵ_d 能级杂质线性耦合到连续谱费米浴，U = 0）
#       浴谱 f(D, ϵ) = sqrt(1-(ϵ/D)^2)/π，[-D, D]，与 test/util.jl 一致
#
# 参考解：ImpurityModelBase 解析解
#   - 虚时间 G(τ) = <d(τ) d†(0)>_β：toulouse_Gτ（相互作用热态）
#   - 实时间 G>(t) = -i<d(t) d†(0)>：toulouse_Gt（直积初态版本）
#
# 运行方式：
#   OMP_NUM_THREADS=1 julia --project=. performance/TTIIF/imaginarytime.jl
#
# 注意：影响泛函构建计时的对照基准不同 ——
#   PartialIF 的基准是 2^k 次 partial-MPO 乘法本身（hybriddynamics 即为全过程），
#   XTRGIF / TDVPIF 同理，hybriddynamics 就包含其全部工作（乘法/求和压缩/流）。
# ============================================================================

using GTEMPO, ImpurityModelBase, LinearAlgebra, Printf, Random

# ------------------------------ 模型参数 -----------------------------------
N = 25
δτ = 0.01
β = N * δτ
ϵ_d = 1.25 * π
D = 10.0
μ = 0.0

f(D, ϵ) = sqrt(1 - (ϵ / D)^2) / π
spec = spectrum(ϵ -> f(D, ϵ), lb = -D, ub = D)

trunc = truncdimcutoff(D = 100, ϵ = 1.0e-8, add_back = 0)

op = [0.0 0; 1 0]   # d

# 四种待比较的影响泛函算法
algs = [
	"PartialIF" => PartialIF(trunc = trunc),
	"XTRGIF"    => XTRGIF(k = 5, fast = true, algmult = SVDCompression(trunc)),
	"ExactTTIIF" => ExactTTIIF(algmult = SVDCompression(trunc)),
	"TDVPIF"    => TDVPIF(trunc = trunc, δ = 0.1),
]

# ------------------------------ 公共函数 -----------------------------------
_relative_error(num, ref) = norm(num - ref) / norm(ref)
_maxbond(mps) = maximum(bond_dimensions(mps))

println("=" ^ 70)
println("Toulouse 模型虚时间基准: N=$N, δτ=$δτ, β=$β, ϵ_d=", round(ϵ_d; digits = 3), ", μ=$μ")
println("trunc: D=$(trunc.D), ϵ=$(trunc.ϵ)")
println("=" ^ 70)

# 解析参考（热态平衡 Green 函数）
bath = fermionicbath(spec, β = β, μ = μ)
b2 = discretebath(bath, δw = 0.2)
exactGτ = toulouse_Gτ(Toulouse(b2, ϵ_d = ϵ_d), collect(0:δτ:β))

# ------------------------------ 各算法计时 ----------------------------------
lattice = GrassmannLattice(N = N, δτ = β / N, contour = :imag, ordering = A1Ā1B1B̄1())
corr = correlationfunction(bath, lattice)

results = Dict{String, Any}()
for (name, alg) in algs
	t_if = @elapsed mpsI = hybriddynamics(lattice, corr, alg)
	mpsI = boundarycondition(mpsI, lattice)

	t_obs = @elapsed begin
		mpsK = sysdynamics(lattice, AndersonIM(μ = ϵ_d, U = 0), trunc = trunc)
		g = Gτ(lattice, mpsK, mpsI)
	end

	err = _relative_error(g, exactGτ)
	println(@sprintf("%-10s IF构建: %6.1f s, 观测量: %5.1f s, 最大键维: %3d, 相对误差: %.3e",
		name, t_if, t_obs, _maxbond(mpsI), err))
	results[name] = (mpsI = mpsI, t_if = t_if, t_obs = t_obs, err = err)
end

# 算法两两之间的 IF 距离（归一化）
names = first.(algs)
for i in 1:length(names), j in (i+1):length(names)
	d = distance(results[names[i]].mpsI, results[names[j]].mpsI) / norm(results[names[i]].mpsI)
	println(@sprintf("distance(%s, %s) = %.3e", names[i], names[j], d))
end
