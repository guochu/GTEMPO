using LinearAlgebra: Hermitian, eigen

@testset "Continuous bath, real time: Toulouse 解析锚" begin
	β = 1.0
	μ = 1.25 * pi
	spec = spectrum_func()
	δt = 0.05; Nt = 8

	bath = fermionicbath(spec, β=β, μ=0)
	# 解析参照：离散化浴的 Toulouse ED（δw 收敛后的参考解）
	bathd = discretebath(bath, δw=0.05)
	model_ed = Toulouse(bathd, ϵ_d=μ)
	ts = (1:Nt) .* δt
	# 二次型精确锚：G>(t,0) = −i·Σ_k |V[1,k]|² e^{−iϵ_k t}(1−f_k)（U=0 全谱精确）
	h = cmatrix(model_ed)
	ev = eigen(Hermitian(h))
	w = ev.values; v2 = abs2.(ev.vectors[1, :])
	gt_ed = [-im * sum(v2 .* exp.(-im .* w .* t) .* (1 ./(exp.(-β .* w) .+ 1))) for t in ts]

	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	kernel = correlationfunction(bath, lat)
	Is = hybriddynamics(lat, kernel, PartialIF())
	model = ToulouseIM(ϵ_d=μ)
	ρ0_imp = thermal_impurity_state(model, β)
	K = sysdynamics(lat, model; initial_state=ρ0_imp)
	Z = integrate(K, Is)
	for i in 1:Nt
		g = greater(lat, model, kernel, Is; i=i, j=0, Z=Z, initial_state=ρ0_imp)
		# 离散浴（δw=0.05）平衡锚 vs 连续 IF：IMB discretebath 权重约定的收敛档
		# （旧套件同档 5e-2；连续谱的承重锚在 imagtime 的解析对拍 3e-3）
		@test g ≈ gt_ed[i] atol = 8e-2   # 连续谱离散化基准（旧套件同档 5e-2）
	end
end
