# kernel 提取 vs QuAPI 直读的交叉核对（IF 后端切换后升级为常驻位级闸门——
# IFD 核 vs QuAPI 原块逐位 ==，后端切换中性证明 + IFD 演化回归锚）。
# QuAPI 侧显式限定（QuAPI 已从 deps 降级 extras——非限定名现绑定 IFD，
# 不限定即退化为自比对恒真，F6 陷阱）。混合 τ↔rt 构造期折叠断言保留。
import QuAPI
@testset "QuAPI 交叉核对：kernel blocks = QuAPI 原块" begin
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	k(b1, b2) = (ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 1, false))

	# 实时：四块直读（canonical = raw；跨支 Keldysh 号 = 指数系数，电流端公式层处理）
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = correlationfunction(bath, lat)
	qt = QuAPI.Δt(bath; N=lat.Nt, t=lat.Nt * lat.δt)   # t 与 correlationfunction 同源（浮点位型一致）
	for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
		M = QuAPI.branch(qt; b1=b1, b2=b2)
		@test kernel.blocks[k(b1, b2)] == M[1:4, 1:4]
	end

	# 虚时：Δτ 块
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = correlationfunction(bath, latτ)
	qτ = QuAPI.Δτ(bath; N=4, δτ=0.25)
	@test kernelτ.blocks[k(:τ, :τ)] == qτ.data[1:4, 1:4]

	# 混合：九块 + τ↔rt 交叉块内嵌负号构造期折叠（旧 _NegView 预翻的物化）
	latm = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernelm = correlationfunction(bath, latm)
	# t/δ 与 correlationfunction 同源（注意 3*0.1 ≠ 0.3——不同位型两次 quadgk 装配差 ~1e-17）
	qm = QuAPI.Δm(bath; Nτ=latm.Nτ, t=latm.Nt * latm.δt, Nt=latm.Nt, δτ=latm.δτ)
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		M = QuAPI.branch(qm; b1=b1, b2=b2)
		r1 = b1 === :τ ? (1:4) : (1:3)
		r2 = b2 === :τ ? (1:4) : (1:3)
		expect = ((b1 === :τ) !== (b2 === :τ)) ? -M[r1, r2] : M[r1, r2]
		@test kernelm.blocks[k(b1, b2)] == expect
	end
end
