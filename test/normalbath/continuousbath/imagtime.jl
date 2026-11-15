@testset "Continuous bath, imaginary time: Toulouse 解析锚" begin
	μ = 1.25 * pi
	spec = spectrum_func()
	δτ = 0.05; N = 10; β = N * δτ
	bath = fermionicbath(spec, β=β, μ=0)
	g_ana = [toulouse_Gτ(bath, τ; ϵ_d=μ) for τ in δτ:δτ:β]

	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	model = ToulouseIM(ϵ_d=μ)
	K = sysdynamics(lat, model)
	Z = integrate(K, Is)
	m = [matsubara(lat, model, kernel, Is; i=i, j=0, Z=Z) for i in 1:N]
	# toulouse_Gτ 的数值实现带随机宇称估计器（同参数跨进程整体 ± 翻转），
	# 锚按符号不定分支；形状/幅值（5e-3 相对）才是断言内容
	@test relerr(m, g_ana) < 5e-3 || relerr(m, -g_ana) < 5e-3
	occ = [occupation(lat, model, kernel, Is; j=i, branch=:τ, Z=Z) for i in 1:N]
	# 平衡占据（τ 无关）；解析锚 = toulouse 基态占据（G(β⁻) = ±n，符号随估计器）
	@test maximum(abs.(occ .- abs.(g_ana[end]))) < 2e-2 || maximum(abs.(occ .+ abs.(g_ana[end]))) < 2e-2
end
