# bridge 层（GMPS 皮）结构查询与算术面：旧 grassmannmps/grassmanntensor/mpo/mult
# 四文件的 Z2Tensors 栈测试已随栈退役；GM 引擎面由 GrassmannMPS 自有套件覆盖。
@testset "API: bridge 结构查询与算术" begin
	lat = GrassmannLattice(N=3, δτ=0.1, bands=2, contour=:imag, symmetry=:u1u1)
	spec = symmetry_spec(:u1u1)
	spaces = site_spaces(spec, length(lat) ÷ 2, 2)
	x = randomgmps(ComplexF64, spaces; D=4)
	@test check_charge(x)
	@test length(bonddims(x)) == length(x) - 1
	# 荷维度：键荷标签 = U1U1 元组
	cd = charge_dimensions(x)
	@test all(d -> all(k -> k isa Tuple{Int,Int}, keys(d)), cd)

	# 缩放：logscale/setlogscale!（GM 直用；scaling/setscaling! 皮肤已删）
	setlogscale!(x, log(2.5))
	@test exp(logscale(x)) ≈ 2.5 atol = 1e-12

	# reorder! 恒等置换 = 值不变（截断关闭）
	y = deepcopy(x)
	reorder!(y, collect(1:length(y)))
	@test scaled_ratio(integrate_scaled(y), integrate_scaled(x)) ≈ 1.0 atol = 1e-12

	# swap! 相邻两站（费米换序号保真：值以链级标度不变量对拍）
	z = deepcopy(x)
	swap!(z, 2)
	@test scaled_ratio(integrate_scaled(z), integrate_scaled(x)) ≈ 1.0 atol = 1e-10

	# multiply（带截断）与朴素楔积一致到截断水平
	w = randomgmps(ComplexF64, spaces; D=2)
	exact = x * w
	truncd = multiply(x, w; trunc=Truncation(maxdim=16, rtol=1e-12))
	@test scaled_ratio(integrate_scaled(truncd), integrate_scaled(exact)) ≈ 1.0 atol = 1e-9
end
