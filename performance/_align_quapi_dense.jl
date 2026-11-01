# 快速对齐测试（N=3）：PartialIF / XTRGIF / ExactTTIIF / TDVPIF 全部与
# QuAPI.Dense 的稠密 IF 逐元对照（_gt_vs_quapi 复用 test/api/quapi_crosscheck.jl）
using Test, Random
using GTEMPO
using Z2Tensors
using ImpurityModelBase
using QuAPI
using QuAPI.Dense

Random.seed!(12354)
include(joinpath(@__DIR__, "..", "test", "api", "quapi_crosscheck.jl"))

@testset "快速对齐: 全算法 vs QuAPI.Dense (N=3)" begin
	β, N = 0.3, 3
	trunc = truncdimcutoff(D=60, ϵ=1.0e-12)
	bath = fermionicbath(DiracDelta(ω=1, α=0.5), β=β)
	tol = 5.0e-2   # TTI 算法含 Prony 截断；PartialIF 应 ~1e-12

	# 虚时
	lat = GrassmannLattice(N=N, δτ=β/N, contour=:imag)
	corr = correlationfunction(bath, lat)
	for (name, alg) in (("PartialIF", PartialIF(trunc=trunc)),
						("XTRGIF", XTRGIF(k=2, alg_mult=SVDCompression(trunc=trunc), verbosity=0)),
						("ExactTTIIF", ExactTTIIF(alg_mult=SVDCompression(trunc=trunc), verbosity=0)),
						("TDVPIF", TDVPIF(trunc=trunc, δ=0.1, verbosity=0)))
		mps = hybriddynamics(lat, corr, alg)
		e = _gt_vs_quapi(_gmps_to_dense(mps), lat, corr)
		@test e < tol
		println("   imag ", name, ": relerr=", e)
	end
	# 实时
	latr = GrassmannLattice(N=N, δt=0.1, contour=:real)
	corrr = correlationfunction(bath, latr)
	for (name, alg) in (("PartialIF", PartialIF(trunc=trunc)),
						("XTRGIF", XTRGIF(k=2, alg_mult=SVDCompression(trunc=trunc), verbosity=0)),
						("ExactTTIIF", ExactTTIIF(alg_mult=SVDCompression(trunc=trunc), verbosity=0)),
						("TDVPIF", TDVPIF(trunc=trunc, δ=0.1, verbosity=0)))
		mps = hybriddynamics(latr, corrr, alg)
		e = _gt_vs_quapi(_gmps_to_dense(mps), latr, corrr)
		@test e < tol
		println("   real ", name, ": relerr=", e)
	end
end
