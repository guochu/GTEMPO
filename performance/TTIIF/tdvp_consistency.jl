# ============================================================================
# 小 N 无截断一致性检验：TDVPIF vs PartialIF（Toulouse 模型，虚时间）
#
# 取 N 足够小使得键维达到可行性上限 2^min(i, L-i) ≪ D（无任何截断），此时
# TDVP 的流形投影恒等，两种算法都给出同一离散化下的精确路径积分，结果应当
# 严格一致（只差 Krylov 求解容差）。
#
# 运行方式：
#   OMP_NUM_THREADS=1 julia --project=. performance/TTIIF/tdvp_consistency.jl
# ============================================================================

using GTEMPO, ImpurityModelBase, LinearAlgebra, Printf, Random

# ------------------------------ 模型参数 -----------------------------------
N = 4              # L = 2N = 8 个格点，最大可行键维 2^4 = 16
δτ = 0.01
β = N * δτ
ϵ_d = 1.25 * π
D = 10.0
μ = 0.0

f(D, ϵ) = sqrt(1 - (ϵ / D)^2) / π
spec = spectrum(ϵ -> f(D, ϵ), lb = -D, ub = D)

# 截断阈值取到机器精度水平、键维上限远大于可行性上限：保证全程无截断
trunc = truncdimcutoff(D = 1000, ϵ = 1.0e-16, add_back = 0)

algs = [
	"PartialIF" => PartialIF(trunc = trunc),
	"TDVPIF"    => TDVPIF(trunc = trunc, δ = 0.1),
]

_relative_error(num, ref) = norm(num - ref) / norm(ref)
_maxbond(mps) = maximum(bonddims(mps))

println("=" ^ 70)
println("无截断一致性检验: N=$N, δτ=$δτ, β=$β (L=$(2N) sites, 可行键维上限 16)")
println("=" ^ 70)

bath = fermionicbath(spec, β = β, μ = 0.0)
b2 = discretebath(bath, δw = 0.2)
exactGτ = toulouse_Gτ(Toulouse(b2, ϵ_d = ϵ_d), collect(0:δτ:β))

lattice = GrassmannLattice(N = N, δτ = β / N, contour = :imag, ordering = A1Ā1B1B̄1())
corr = correlationfunction(bath, lattice)

results = Dict{String, Any}()
for (name, alg) in algs
	t_if = @elapsed mpsI = hybriddynamics(lattice, corr, alg)
	mpsI = boundarycondition(mpsI, lattice)

	t_obs = @elapsed begin
		mpsK = sysdynamics(lattice, AndersonIM(ϵ_d = ϵ_d, U = 0), trunc = trunc)
		g = cached_gf_fast(lattice, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ)
	end

	err = _relative_error(g, exactGτ)
	println(@sprintf("%-10s IF构建: %6.1f s, 观测量: %5.1f s, 最大键维: %3d, GF相对误差: %.3e",
		name, t_if, t_obs, _maxbond(mpsI), err))
	results[name] = (mpsI = mpsI, g = g, err = err)
end

d_if = distance(results["PartialIF"].mpsI, results["TDVPIF"].mpsI) / norm(results["PartialIF"].mpsI)
d_g = _relative_error(results["TDVPIF"].g, results["PartialIF"].g)
println(@sprintf("IF  相对距离 (PartialIF vs TDVPIF): %.3e", d_if))
println(@sprintf("GF  相对偏差 (PartialIF vs TDVPIF): %.3e", d_g))
