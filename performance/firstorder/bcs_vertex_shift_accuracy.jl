# BCS bath 的顶点位移 A/B 精度基准 —— 与 vertex_shift_accuracy.jl（普通 bath）
# 配套：BCS 访问器目前 τ 环/RT/混合全部 raw（"BCS RT 暂不迁移"），本基准
# 回答：BCS 是否也呈现普通 bath 的三轴格局（τ 环 vertex 劣 ~30×、实时 RT
# vertex 优 ~25×）？
#
# - raw：当前语义（BCS getindex 三分支全部原始读取）；
# - vertex：τ 环与实时 RT 施加与普通 bath 相同的顶点位移（_vs_τ_acc /
#   _vs_rt_acc），混合轴维持 raw（2026-10-30 决议：混合不位移）。
#
# 参考：few-mode BCS ED（bcs_ed + gτ_ed / greater_lesser_ed，与
# test/bcs/fewmodesbath 同源）。
#
# 运行：julia --project=../.. performance/firstorder/bcs_vertex_shift_accuracy.jl
#
# 实测结论（2026-10-30，AndersonIM U=1, ϵ_d=−0.8, Δ=0.6, DiracDelta(ω=1, α=0.5),
# β=1，PartialIF，D=60；参考 few-mode bcs_ed）：
#
#   虚时 τ 环（Nτ=8..32）：vertex 0.017→0.0076 (p≈0.6)
#                         raw    0.0049→0.0046（平底 —— few-mode BCS ED 参考
#                                自身带 ~0.48% 残差：D=60 与 D=250 逐位不变，
#                                连续谱参考下同管线仅 8e-4）
#                         → vertex 差 1.7-3.5×：与普通 bath 同向，τ 环维持 raw ✓
#   实时 RT（Nt=4..16）：vertex gt 2.2e-4 平底（截断地板）；raw 0.0036→0.0008
#                        (p≈1.1) → vertex 好 4-16×：与普通 bath 的 25× 同向。
#                        决议（2026-10-30）：BCS RT 已迁移 vertex（迁移后
#                        gt/lt 实测 4.2e-4/2.1e-4，raw 约 4.2e-3/2.4e-3）；
#                        QuAPI.Dense 的 BCS RT 仍为旧装配，待同步。
#   另：τ 环回退 raw 后，BCS continuous 复能隙系统残差 0.0254（δw 不敏感）
#   → 8.2e-4 —— 此前归因"反常块约定待定"的残差实为 τ 环位移所致，已消。

using Random, LinearAlgebra
using GTEMPO
using Z2Tensors
using ImpurityModelBase
using QuAPI
using QuAPI: GrassmannContourIndex, branch, branches, _vs_τ_acc, _vs_rt_acc

Random.seed!(12354)

# ED 参考构造器来自测试公用工具
using GTEMPO: A1B1B̄1Ā1, Ā2A1B̄2B1,
          A1B1ā1b̄1Ā1B̄1a1b1,
          Ā2A1ā1a2B̄2B1b̄1b̄2,
          Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2, A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1,
          FockMatrix, fock_propagator, fock_thermalstate
const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "test", "util.jl"))
include(joinpath(GT_ROOT, "test", "normalbath", "util.jl"))

# ---------------------------------------------------------------------------
# 两种 BCS 访问器语义（顶层逐条切换；勿在函数内安装——见主基准脚本的
# world-age 注记）
# ---------------------------------------------------------------------------
function install_bcs_raw!()
	@eval QuAPI function Base.getindex(x::BCSCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		σᵢ = isodd(i.band) ? :↑ : :↓
		σⱼ = isodd(j.band) ? :↑ : :↓
		Δb = (σᵢ === σⱼ) ? x.uu : (σᵢ === :↑ ? x.ud : x.du)
		v = branch(Δb; b1=branch(i), b2=branch(j))[i.j, j.j]
		return (i.conj && !j.conj) ? -v : v
	end
end

# 2026-10-30 起的线上语义：τ 环 raw、实时 RT vertex、混合 raw
function install_bcs_live!()
	@eval QuAPI function Base.getindex(x::BCSCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		σᵢ = isodd(i.band) ? :↑ : :↓
		σⱼ = isodd(j.band) ? :↑ : :↓
		Δb = (σᵢ === σⱼ) ? x.uu : (σᵢ === :↑ ? x.ud : x.du)
		z = zero(eltype(branch(Δb; b1=branch(i), b2=branch(j))))
		bs = branches(x)
		if bs == (:τ,)
			v = branch(Δb; b1=:τ, b2=:τ)[i.j, j.j]
		elseif bs == (:+, :-)
			m = _vs_rt_acc(branch(i), i.j, branch(j), j.j,
						   length(branch(x.uu; b1=:+, b2=:+).ηⱼₖ) - 1)
			m === nothing && return z
			v = m[5] * branch(Δb; b1=m[1], b2=m[3])[m[2], m[4]]
		else
			v = branch(Δb; b1=branch(i), b2=branch(j))[i.j, j.j]
		end
		return (i.conj && !j.conj) ? -v : v
	end
end

function install_bcs_vertex!()
	@eval QuAPI function Base.getindex(x::BCSCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		σᵢ = isodd(i.band) ? :↑ : :↓
		σⱼ = isodd(j.band) ? :↑ : :↓
		Δb = (σᵢ === σⱼ) ? x.uu : (σᵢ === :↑ ? x.ud : x.du)
		bs = branches(x)
		z = zero(eltype(branch(Δb; b1=branch(i), b2=branch(j))))
		if bs == (:τ,)
			Nτ = length(branch(x.uu; b1=:τ, b2=:τ).ηⱼₖ)
			m = Nτ == 1 ? (:τ, 1, :τ, 1, 1) : _vs_τ_acc(i.j, j.j, Nτ)
			v = m[5] * branch(Δb; b1=m[1], b2=m[3])[m[2], m[4]]
		elseif bs == (:+, :-)
			m = _vs_rt_acc(branch(i), i.j, branch(j), j.j,
						   length(branch(x.uu; b1=:+, b2=:+).ηⱼₖ) - 1)
			m === nothing && return z
			v = m[5] * branch(Δb; b1=m[1], b2=m[3])[m[2], m[4]]
		else
			# 混合维持 raw（2026-10-30 决议），本 A/B 不涉及
			v = branch(Δb; b1=branch(i), b2=branch(j))[i.j, j.j]
		end
		return (i.conj && !j.conj) ? -v : v
	end
end

# ---------------------------------------------------------------------------
# 测量（few-mode BCS ED 参考）
# ---------------------------------------------------------------------------
const Ub, ϵb, ωb, αb, Δb_gap = 1.0, 0.8, 1.0, 0.5, 0.6
const βb = 1.0

function bcs_imag_error(Nτ::Int; trunc)
	δτ = βb / Nτ
	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
	model = AndersonIM(U=Ub, ϵ_d=-ϵb)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	corr = correlationfunction(bcsbath(fermionicbath(DiracDelta(ω=ωb, α=αb), β=βb), Δ=Δb_gap), lattice)
	mpsI = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
	g = gτ_series(lattice, mpsK, mpsI)
	H, a, adag, H0 = bcs_ed(Ub, ϵb, ω₀=ωb, α=αb, Δ=Δb_gap)
	ref = gτ_ed(H, a, adag, collect(0:δτ:βb), βb)
	return norm(g - ref) / norm(ref)
end

function bcs_real_error(Nt::Int; trunc, T=0.2)
	δt = T / Nt
	ts = collect(0:δt:T)
	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=2)
	model = AndersonIM(U=Ub, ϵ_d=-ϵb)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=βb)
	corr = correlationfunction(bcsbath(fermionicbath(DiracDelta(ω=ωb, α=αb), β=βb), Δ=Δb_gap), lattice)
	mpsI = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI)
	g1 = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	g2 = [-im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	H, a, adag, H0 = bcs_ed(Ub, ϵb, ω₀=ωb, α=αb, Δ=Δb_gap)
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, βb)
	return norm(g1 - gt_ed) / norm(gt_ed), norm(g2 - lt_ed) / norm(lt_ed)
end

_m(v) = sum(v) / length(v)
slope(hs, es) = (sum((hs .- _m(hs)) .* log.(es)) / sum(abs2, hs .- _m(hs))) / log(2)

# ---------------------------------------------------------------------------
# 顶层批次（world-age：见主基准脚本注记）
# ---------------------------------------------------------------------------
trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

println("=" ^ 78)
println("BCS 虚时 τ 环：g(τ) vs bcs_ed   (β=$βb, Δ=$Δb_gap, AndersonIM U=$Ub, D=60)")
println("=" ^ 78)
Ns = [8, 16, 32]
println(rpad("Nτ", 6), rpad("δτ", 10), rpad("vertex relerr", 16), rpad("raw relerr", 16), "vertex/raw")
install_bcs_vertex!()
bcs_imag_vertex = [bcs_imag_error(N; trunc=trunc) for N in Ns]
install_bcs_raw!()
bcs_imag_raw = [bcs_imag_error(N; trunc=trunc) for N in Ns]
for (k, N) in enumerate(Ns)
	println(rpad(N, 6), rpad(round(βb / N; sigdigits=3), 10),
			rpad(round(bcs_imag_vertex[k]; sigdigits=3), 16),
			rpad(round(bcs_imag_raw[k]; sigdigits=3), 16),
			round(bcs_imag_vertex[k] / bcs_imag_raw[k]; sigdigits=3))
end
hs = -log2.(Ns)
println("收敛阶：vertex p = ", round(slope(hs, bcs_imag_vertex); digits=2),
		"    raw p = ", round(slope(hs, bcs_imag_raw); digits=2))
println()

println("=" ^ 78)
println("BCS 实时 Keldysh：gt/lt vs bcs_ed   (T=0.2, β=$βb, Δ=$Δb_gap, D=60)")
println("=" ^ 78)
Ns = [4, 8, 16]
println(rpad("Nt", 6), rpad("δt", 10), rpad("vertex gt|lt", 24), rpad("raw gt|lt", 24), "vertex/raw (gt)")
install_bcs_vertex!()
bcs_real_vertex = [bcs_real_error(N; trunc=trunc) for N in Ns]
install_bcs_raw!()
bcs_real_raw = [bcs_real_error(N; trunc=trunc) for N in Ns]
for (k, N) in enumerate(Ns)
	v, r = bcs_real_vertex[k][1], bcs_real_raw[k][1]
	vl, rl = bcs_real_vertex[k][2], bcs_real_raw[k][2]
	println(rpad(N, 6), rpad(round(0.2 / N; sigdigits=3), 10),
			rpad("$(round(v; sigdigits=4)) | $(round(vl; sigdigits=4))", 24),
			rpad("$(round(r; sigdigits=4)) | $(round(rl; sigdigits=4))", 24),
			round(v / r; sigdigits=3))
end
hs = -log2.(Ns)
println("gt 收敛阶：vertex p = ", round(slope(hs, [x[1] for x in bcs_real_vertex]); digits=2),
		"    raw p = ", round(slope(hs, [x[1] for x in bcs_real_raw]); digits=2))
println()

# ---- 存盘：CSV + 两面板 SVG（与主基准同款手写绘图）----
const OUTDIR = @__DIR__
open(joinpath(OUTDIR, "bcs_vertex_shift_error_scaling.csv"), "w") do io
	println(io, "axis,N,delta,convention,relerr")
	for (k, N) in enumerate([8, 16, 32])
		println(io, "BCS imaginary,N,", βb / N, ",vertex,", bcs_imag_vertex[k])
		println(io, "BCS imaginary,N,", βb / N, ",raw,", bcs_imag_raw[k])
	end
	for (k, N) in enumerate([4, 8, 16])
		println(io, "BCS real,N,", 0.2 / N, ",vertex,", bcs_real_vertex[k][1])
		println(io, "BCS real,N,", 0.2 / N, ",raw,", bcs_real_raw[k][1])
	end
end

function svg_panel(x0, y0, W, H, title, xs, ev, er, xlab, note)
	xmin, xmax = floor(log10(minimum(xs))), ceil(log10(maximum(xs)))
	allof = vcat(ev, er); ymin, ymax = floor(log10(minimum(allof))), ceil(log10(maximum(allof)))
	ml, mr, mt, mb = 54, 12, 26, 40
	px(x) = x0 + ml + (log10(x) - xmin) / (xmax - xmin) * (W - ml - mr)
	py(y) = y0 + H - mb - (log10(y) - ymin) / (ymax - ymin) * (H - mt - mb)
	poly(es) = join(["$(round(px(xs[k]); digits=1)),$(round(py(es[k]); digits=1))" for k in eachindex(es)], " ")
	io = IOBuffer()
	write(io, """<rect x="$(x0)" y="$(y0)" width="$(W)" height="$(H)" fill="white" stroke="#999"/>
<text x="$(x0+W/2)" y="$(y0+15)" text-anchor="middle" font-size="13" font-family="serif" font-style="italic">$title</text>
<text x="$(x0+W/2)" y="$(y0+H-6)" text-anchor="middle" font-size="11" font-family="serif">$xlab</text>
<text x="$(x0+W/2)" y="$(y0+28)" text-anchor="middle" font-size="9" fill="#555">$note</text>
""")
	for t in Int(xmin):Int(xmax)
		xx = round(px(10.0^t); digits=1)
		write(io, "<line x1='$xx' y1='$(y0+mt)' x2='$xx' y2='$(y0+H-mb)' stroke='#e8e8e8'/><text x='$xx' y='$(y0+H-mb+13)' text-anchor='middle' font-size='9'>10<tspan dy='-4' font-size='7'>$t</tspan></text>\n")
	end
	for t in Int(ymin):Int(ymax)
		yy = round(py(10.0^t); digits=1)
		write(io, "<line x1='$(x0+ml)' y1='$yy' x2='$(x0+W-mr)' y2='$yy' stroke='#e8e8e8'/><text x='$(x0+ml-5)' y='$(yy+3)' text-anchor='end' font-size='9'>10<tspan dy='-4' font-size='7'>$t</tspan></text>\n")
	end
	write(io, "<polyline points='$(poly(ev))' fill='none' stroke='#d62728' stroke-width='1.6'/>
")
	write(io, "<polyline points='$(poly(er))' fill='none' stroke='#1f77b4' stroke-width='1.6' stroke-dasharray='7 3'/>
")
	for k in eachindex(xs)
		write(io, "<circle cx='$(round(px(xs[k]);digits=1))' cy='$(round(py(ev[k]);digits=1))' r='3' fill='#d62728'/><circle cx='$(round(px(xs[k]);digits=1))' cy='$(round(py(er[k]);digits=1))' r='3' fill='#1f77b4'/>
")
	end
	String(take!(io))
end
W, H, GAP = 360, 300, 18
svg = ["<svg xmlns='http://www.w3.org/2000/svg' width='$(2W+GAP+20)' height='$(H+52)' font-family='serif'>",
	   "<rect width='100%' height='100%' fill='white'/>",
	   "<text x='$(W+GAP/2+10)' y='20' text-anchor='middle' font-size='15'>BCS IF vertex-shift vs raw: error scaling (PartialIF, bcs_ed)</text>",
	   svg_panel(10, 32, W, H, "BCS imaginary (τ ring)", [βb/N for N in [8,16,32]], bcs_imag_vertex, bcs_imag_raw, "δτ", "raw floor ~5e-3 = ED-reference residual (D-independent)"),
	   svg_panel(10 + W + GAP, 32, W, H, "BCS real (Keldysh, gt)", [0.2/N for N in [4,8,16]], [x[1] for x in bcs_real_vertex], [x[1] for x in bcs_real_raw], "δt", "vertex floor ~2e-4 = truncation floor"),
	   "<line x1='$(W-80)' y1='$(H+42)' x2='$(W-60)' y2='$(H+42)' stroke='#d62728' stroke-width='1.6'/><text x='$(W-54)' y='$(H+46)' font-size='11'>vertex</text>",
	   "<line x1='$(W+20)' y1='$(H+42)' x2='$(W+40)' y2='$(H+42)' stroke='#1f77b4' stroke-width='1.6' stroke-dasharray='7 3'/><text x='$(W+46)' y='$(H+46)' font-size='11'>raw (cell)</text>",
	   "</svg>"]
open(joinpath(OUTDIR, "bcs_vertex_shift_error_scaling.svg"), "w") do io
	write(io, join(svg, "\n"))
end
println("已保存：bcs_vertex_shift_error_scaling.svg / .csv（performance/firstorder/）")

install_bcs_live!()   # 恢复当前语义
println("（BCS 访问器已恢复为当前语义：τ 环/混合 raw、实时 RT vertex）")
