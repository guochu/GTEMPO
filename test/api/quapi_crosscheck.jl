# kernel 提取 vs QuAPI 直读的交叉核对（T7 退役后唯一保留的 QuAPI 面）
# （QuAPI 经 @reexport 直接可用）
@testset "QuAPI 交叉核对：kernel blocks = QuAPI 原块" begin
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)

	# 实时：四块直读（canonical = raw；跨支 Keldysh 号由消费端处理）
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = hybridization(lat, bath)
	qt = Δt(bath; N=4, t=0.4)
	for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
		M = branch(qt; b1=b1, b2=b2)
		@test kernel.blocks[(b1, b2)] ≈ M[1:4, 1:4] atol = 1e-14
	end

	# 虚时：Δτ 块
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = hybridization(latτ, bath)
	qτ = Δτ(bath; N=4, δτ=0.25)
	@test kernelτ.blocks[(:τ, :τ)] ≈ qτ.data[1:4, 1:4] atol = 1e-14

	# 混合：九块 + τ↔rt 交叉块内嵌负号预翻（KernelCorrView 语义）
	latm = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernelm = hybridization(latm, bath)
	qm = Δm(bath; Nτ=4, t=0.3, Nt=3, δτ=0.25)
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		M = branch(qm; b1=b1, b2=b2)
		r1 = b1 === :τ ? (1:4) : (1:3)
		r2 = b2 === :τ ? (1:4) : (1:3)
		expect = ((b1 === :τ) !== (b2 === :τ)) ? -M[r1, r2] : M[r1, r2]
		@test kernelm.blocks[(b1, b2)] ≈ expect atol = 1e-14
	end
end
