@testset "API: correlationfunction（QuAPI 包装）" begin
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)

	# 虚时：Δτ 矩阵
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	corr = correlationfunction(bath, lat)
	@test corr isa ImagCorrelationFunction
	# 实时：QuAPI 关联对象四块
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	corr = correlationfunction(bath, lat)
	@test corr isa RealCorrelationFunction
	# 混合
	lat = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	corr = correlationfunction(bath, lat)
	@test corr isa MixedCorrelationFunction
end

@testset "API: baths（ImpurityModelBase 皮肤）" begin
	β = 1.0
	bath = discretefermionicbath([1.0, 2.0], [0.5, 0.3]; β=β)
	@test bath.β == β
	bath = fermionicbath(spectrum_func(), β=β, μ=0.3)
	@test bath.β == β && bath.μ == 0.3
	@test bath isa AbstractNormalBath
end
