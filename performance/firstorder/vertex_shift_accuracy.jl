# 顶点位移装配 vs 原始（cell）读取的 GTEMPO 精度基准 —— Toulouse 模型
#
# 回答的问题：QuAPI 2026-10-27 起在关联函数访问器里引入的「顶点位置装配」
# （τ 环 skew-circulant 重组、实时同支块列移 + 转折点开窗），相对之前的
# cell 直读约定，在 GTEMPO 的影响泛函精度上到底有没有优势？
#
# 方法：
# - vertex：把 QuAPI 当前访问器语义（contourindex.jl 的
#   getindex(::AbstractCorrelationFunction, ::GrassmannContourIndex, ...)，
#   虚时/实时分支）原样复制进本脚本，可随时切换回去；
# - raw：branch(x; b1, b2)[i.j, j.j] 直读（2026-10-27 之前的 cell 约定）。
#   实时跨支块的 Keldysh 内嵌符号由 GTEMPO 的 _physical_coef 统一补偿，
#   两种约定下补偿语义一致，因此 A/B 的差异纯粹是「位移 vs 不位移」。
# - 两轴分别测量（顶点位移在这两轴都处于激活状态；mixed 已于 2026-10-29
#   回退 raw，不在此基准内）：
#   * 虚时 τ 环：g(τ) vs 解析 toulouse_Gτ，β 固定，Nτ ∈ {8,16,32,64}；
#   * 实时 Keldysh：gt(t) vs 解析 im·toulouse_Gt，T 固定，Nt ∈ {4,8,16,32}。
# - 截断取得足够紧（D=200, ϵ=1e-12），使测到的误差是离散化误差而非键截断。
#
# 运行：julia --project=../.. performance/firstorder/vertex_shift_accuracy.jl

# ---------------------------------------------------------------------------
# 实测结论（2026-10-30，DiracDelta(ω=1, α=0.5), β=1, ϵ_d=0.7, T=0.4,
# PartialIF 静态构造，D=200/ϵ=1e-12，单模 ED 参考）：
#
#   虚时 τ 环（g(τ),  Nτ=8..64）：vertex 0.017→0.0021 (p≈1.0)
#                                raw    0.0006→5e-5 (p≈1.1)
#                                → 顶点位移差 ~30 倍（同为低阶收敛，常数劣势）
#   实时 Keldysh（gt/lt, Nt=4..32）：vertex 0.00052→6e-5 (p≈1.0)
#                                raw    0.014→0.0015 (p≈1.1)
#                                → 顶点位移好 ~25 倍
#   （混合轴的顶点位移已被 2026-10-29 的 ED 判决回退：6-7% 且不随 δt 收敛）
#
#   混合 Kadanoff（gt, N=3..24）：vertex 0.079→0.056 **不收敛**（p≈0.17，
#                                卡在 5-8% 系统偏差——2026-10-29 ED 判决的量化）
#                                raw    0.009→0.0012 (p≈0.97)，gτ 达 p≈1.4
#
# 决议（2026-10-30，用户拍板）：实时 Keldysh 保留顶点位移；τ 环与混合轴
# 回退 raw（QuAPI 访问器与 GTEMPO TTI :τ bulk 已同步回退，全套测试绿）。
# QuAPI.Dense 对齐测试继续屏蔽，以严格解（ED/解析）为准确认约定。
# ---------------------------------------------------------------------------

using Random, LinearAlgebra
using GTEMPO
using Z2Tensors
using ImpurityModelBase
using QuAPI
using QuAPI: GrassmannContourIndex, branch, branches, _vs_τ_acc, _vs_rt_acc

# ED 参考构造器来自测试公用工具（greater_lesser_ed / singlemode_ed）
using GTEMPO: A1B1B̄1Ā1, Ā2A1B̄2B1,
          A1B1ā1b̄1Ā1B̄1a1b1,
          Ā2A1ā1a2B̄2B1b̄1b̄2,
          Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2, A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1,
          FockMatrix, fock_propagator,
          fock_thermalstate
const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "test", "util.jl"))
include(joinpath(GT_ROOT, "test", "normalbath", "util.jl"))

Random.seed!(12354)

# ---------------------------------------------------------------------------
# 两种访问器约定的安装器（签名与 QuAPI 侧一致，后定义者生效）
# ---------------------------------------------------------------------------
function install_vertex_accessor!()
	@eval QuAPI function Base.getindex(x::AbstractCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		(i.conj == j.conj) && return zero(eltype(branch(x; b1=branch(i), b2=branch(j))))
		bs = branches(x)
		m = if bs == (:τ,)
			Nτ = length(branch(x; b1=:τ, b2=:τ).ηⱼₖ)
			Nτ == 1 ? (:τ, 1, :τ, 1, 1) : _vs_τ_acc(i.j, j.j, Nτ)
		elseif bs == (:+, :-)
			_vs_rt_acc(branch(i), i.j, branch(j), j.j,
					   length(branch(x; b1=:+, b2=:+).ηⱼₖ) - 1)
		else
			error("vertex accessor benchmark covers the imag/real contours only")
		end
		m === nothing && return zero(eltype(branch(x; b1=branch(i), b2=branch(j))))
		D = branch(x; b1=m[1], b2=m[3])
		v = m[5] * D[m[2], m[4]]
		return i.conj ? -v : v
	end
end

function install_raw_accessor!()
	@eval QuAPI function Base.getindex(x::AbstractCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		(i.conj == j.conj) && return zero(eltype(branch(x; b1=branch(i), b2=branch(j))))
		v = branch(x; b1=branch(i), b2=branch(j))[i.j, j.j]
		return i.conj ? -v : v
	end
end

# 混合轴的 A/B：raw 即当前语义（2026-10-29 回退），vertex 为当年被 ED 判决
# 否掉的 `_vs_mixed_acc` 位移版（仅打补丁于 AbstractMixedCorrelationFunction）
#（只依赖 branch() 与分块尺寸，不触内部字段——mixed 结构仍在重构中）
function install_mixed_vertex_accessor!()
	@eval QuAPI function Base.getindex(x::AbstractMixedCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		z = zero(eltype(x.Gττ))
		(i.conj == j.conj) && return z
		m = _vs_mixed_acc(branch(i), i.j, branch(j), j.j,
						  size(x.Gτ₊, 1), size(x.G₊τ, 1) - 1)
		m === nothing && return z
		v = m[5] * branch(x; b1=m[1], b2=m[3])[m[2], m[4]]
		return i.conj ? -v : v
	end
end

function install_mixed_raw_accessor!()
	@eval QuAPI function Base.getindex(x::AbstractMixedCorrelationFunction,
					   i::GrassmannContourIndex, j::GrassmannContourIndex)
		(i.conj == j.conj) && return zero(eltype(x.Gττ))
		v = branch(x; b1=branch(i), b2=branch(j))[i.j, j.j]
		return i.conj ? -v : v
	end
end

# ---------------------------------------------------------------------------
# Toulouse 管线（半圆谱，解析可解）
# ---------------------------------------------------------------------------
const β = 1.0
const ϵ_d = 0.7
const ωb, αb = 1.0, 0.5
toulouse_bath() = fermionicbath(DiracDelta(ω=ωb, α=αb), β=β)

function imag_error(N::Int; trunc)
	# 单模 bath + ED 精确参考（与 fewmodesbath/imagtime.jl 同源）。连续谱的
	# 解析 toulouse_Gτ 标量参考在 ϵ_d 取一般值时其内参约定不可靠（ϵ_d=±0.5
	# 输出逐位相同），不适合做基准参考。
	lat = GrassmannLattice(N=N, δτ=β / N, contour=:imag)
	bath = toulouse_bath()
	corr = correlationfunction(bath, lat)
	H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=0, bathspecs=[(ωb, αb)])
	ref = gτ_ed(H, a, adag, collect(0:β/N:β), β)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, ToulouseIM(ϵ_d=ϵ_d), trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	cache = environments(lat, mpsK, mpsI)
	g = cached_gf_fast(lat, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	return norm(g - ref) / norm(ref)
end

function realtime_error(Nt::Int; trunc, T=0.4)
	# 单模 bath（DiracDelta）+ ED 精确参考（与 fewmodesbath/realtime.jl 同源）；
	# 连续谱的解析 toulouse_Gt 在 ϵ_d 远离大正值时其内参约定不可靠
	#（ϵ_d=±0.5 输出逐位相同），不适合做基准参考。
	δt = T / Nt
	bath = toulouse_bath()
	H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=0, bathspecs=[(ωb, αb)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, collect(0:δt:T), β)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	model = ToulouseIM(ϵ_d=ϵ_d)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
	cache = environments(lat, mpsK, mpsI)
	gt = [-im * cached_greater(lat, k, mpsK, mpsI, cache=cache) for k in 1:lat.k]
	lt = [-im * cached_lesser(lat, k, mpsK, mpsI, cache=cache) for k in 1:lat.k]
	return norm(gt - gt_ed) / norm(gt_ed), norm(lt - lt_ed) / norm(lt_ed)
end

function mixed_error(N::Int; trunc, T=0.2)
	# Kadanoff 轮廓：Nt = Nτ = N 同步加密（δt = T/N, δτ = β/N）；
	# 参考 = 单模 ED（greater_lesser_ed_mixed + gτ_ed，与 fewmodesbath/mixedtime.jl 同源）
	δt, δτ = T / N, β / N
	ts = collect(0:δt:T)
	H, a, adag, H0 = singlemode_ed(μ=ϵ_d, U=0, bathspecs=[(ωb, αb)])
	gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
	gτ_ref = gτ_ed(H, a, adag, collect(0:δτ:β), β)
	lat = GrassmannLattice(Nt=N, δt=δt, Nτ=N, δτ=δτ, contour=:mixed)
	corr = correlationfunction(toulouse_bath(), lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, ToulouseIM(ϵ_d=ϵ_d), trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
	return norm(gt - gt_ed) / norm(gt_ed), norm(gτ - gτ_ref) / norm(gτ_ref)
end

_m(v) = sum(v) / length(v)
# hs 以 log2 为单位、误差取自然对数：除以 log(2) 归一到「误差对 δ 的阶」
slope(hs, es) = (sum((hs .- _m(hs)) .* log.(es)) / sum(abs2, hs .- _m(hs))) / log(2)

# 顶层执行：每条顶层语句都在最新 world 里求值。务必不要把这些批次搬进
# 某个函数里 —— 函数栈帧停留在其启动时的 world，中途替换访问器方法后
# 后续调用仍会用旧方法（这正是本脚本第一版虚/实两轴"结果相同"的原因）。
trunc = truncdimcutoff(D=200, ϵ=1.0e-12)

println("=" ^ 78)
println("虚时 τ 环：g(τ) vs 单模 ED   (β=$β, ϵ_d=$ϵ_d, DiracDelta(ω=$ωb, α=$αb), D_bond=200)")
println("=" ^ 78)
Ns = [8, 16, 32, 64]
println(rpad("Nτ", 6), rpad("δτ", 10), rpad("vertex relerr", 18), rpad("raw relerr", 18), "vertex/raw")
install_vertex_accessor!()
imag_vertex = [imag_error(N; trunc=trunc) for N in Ns]
install_raw_accessor!()
imag_raw = [imag_error(N; trunc=trunc) for N in Ns]
for (k, N) in enumerate(Ns)
	v, r = imag_vertex[k], imag_raw[k]
	println(rpad(N, 6), rpad(round(β / N; sigdigits=3), 10),
			rpad(round(v; sigdigits=3), 18), rpad(round(r; sigdigits=3), 18),
			round(v / r; sigdigits=3))
end
hs = -log2.(Ns)
println("收敛阶：vertex p = ", round(slope(hs, imag_vertex); digits=2),
		"    raw p = ", round(slope(hs, imag_raw); digits=2))
println()

println("=" ^ 78)
println("实时 Keldysh：gt/lt vs 单模 ED   (T=0.4, β=$β, ϵ_d=$ϵ_d, DiracDelta(ω=$ωb, α=$αb))")
println("=" ^ 78)
Ns = [4, 8, 16, 32]
println(rpad("Nt", 6), rpad("δt", 10), rpad("vertex gt|lt", 24), rpad("raw gt|lt", 24), "vertex/raw (gt)")
install_vertex_accessor!()
real_vertex = [realtime_error(N; trunc=trunc) for N in Ns]
install_raw_accessor!()
real_raw = [realtime_error(N; trunc=trunc) for N in Ns]
for (k, N) in enumerate(Ns)
	v, r = real_vertex[k][1], real_raw[k][1]
	vl, rl = real_vertex[k][2], real_raw[k][2]
	println(rpad(N, 6), rpad(round(0.4 / N; sigdigits=3), 10),
			rpad("$(round(v; sigdigits=4)) | $(round(vl; sigdigits=4))", 24),
			rpad("$(round(r; sigdigits=4)) | $(round(rl; sigdigits=4))", 24),
			round(v / r; sigdigits=3))
end
hs = -log2.(Ns)
println("gt 收敛阶：vertex p = ", round(slope(hs, [x[1] for x in real_vertex]); digits=2),
		"    raw p = ", round(slope(hs, [x[1] for x in real_raw]); digits=2))
println("lt 收敛阶：vertex p = ", round(slope(hs, [x[2] for x in real_vertex]); digits=2),
		"    raw p = ", round(slope(hs, [x[2] for x in real_raw]); digits=2))
println()

println("=" ^ 78)
println("混合 Kadanoff：gt/gτ vs 单模 ED   (T=0.2, β=$β, Nt=Nτ=N, DiracDelta(ω=$ωb, α=$αb))")
println("=" ^ 78)
Ns = [3, 6, 12, 24]
println(rpad("N", 6), rpad("δt|δτ", 16), rpad("vertex gt|gτ", 24), rpad("raw gt|gτ", 24), "vertex/raw (gt)")
install_mixed_vertex_accessor!()
mixed_vertex = [mixed_error(N; trunc=trunc) for N in Ns]
install_mixed_raw_accessor!()
mixed_raw = [mixed_error(N; trunc=trunc) for N in Ns]
for (k, N) in enumerate(Ns)
	v, r = mixed_vertex[k][1], mixed_raw[k][1]
	vl, rl = mixed_vertex[k][2], mixed_raw[k][2]
	println(rpad(N, 6), rpad("$(round(0.2/N; sigdigits=3))|$(round(β/N; sigdigits=3))", 16),
			rpad("$(round(v; sigdigits=4)) | $(round(vl; sigdigits=4))", 24),
			rpad("$(round(r; sigdigits=4)) | $(round(rl; sigdigits=4))", 24),
			round(v / r; sigdigits=3))
end
hs = -log2.(Ns)
println("gt 收敛阶：vertex p = ", round(slope(hs, [x[1] for x in mixed_vertex]); digits=2),
		"    raw p = ", round(slope(hs, [x[1] for x in mixed_raw]); digits=2))
println("gτ 收敛阶：vertex p = ", round(slope(hs, [x[2] for x in mixed_vertex]); digits=2),
		"    raw p = ", round(slope(hs, [x[2] for x in mixed_raw]); digits=2))
println()

# ---- 三轴 error scaling 画图（无绘图依赖，手写 SVG）+ 数据存盘 ----
const OUTDIR = @__DIR__
panels = [
	("imaginary (τ ring)", [β / N for N in [8, 16, 32, 64]], imag_vertex, imag_raw, "δτ"),
	("real (Keldysh)", [0.4 / N for N in [4, 8, 16, 32]], [x[1] for x in real_vertex], [x[1] for x in real_raw], "δt"),
	("mixed (Kadanoff, gt)", [0.2 / N for N in [3, 6, 12, 24]], [x[1] for x in mixed_vertex], [x[1] for x in mixed_raw], "δt"),
]
# CSV
const PANEL_NS = [[8, 16, 32, 64], [4, 8, 16, 32], [3, 6, 12, 24]]
open(joinpath(OUTDIR, "vertex_shift_error_scaling.csv"), "w") do io
	println(io, "axis,N,delta,convention,relerr")
	for (k, (name, ds, vs, rs, _)) in enumerate(panels), (tag, es) in (("vertex", vs), ("raw", rs))
		for m in eachindex(ds)
			println(io, name, ",", PANEL_NS[k][m], ",", ds[m], ",", tag, ",", es[m])
		end
	end
end
# SVG: 三个 log-log 面板 + 一阶/二阶参考线
function svg_loglog_panel(x0, y0, W, H, title, xs, ev, er, xlab)
	# 映射（轴界向外对齐到十进制刻度，刻度线就不会越出面板）
	xmin, xmax = floor(log10(minimum(xs))), ceil(log10(maximum(xs)))
	allof = vcat(ev, er)
	ymin, ymax = floor(log10(minimum(allof))), ceil(log10(maximum(allof)))
	ml, mr, mt, mb = 54, 12, 26, 40
	px(x) = x0 + ml + (log10(x) - xmin) / (xmax - xmin) * (W - ml - mr)
	py(y) = y0 + H - mb - (log10(y) - ymin) / (ymax - ymin) * (H - mt - mb)
	poly(es) = join(["$(round(px(xs[k]); digits=1)),$(round(py(es[k]); digits=1))" for k in eachindex(es)], " ")
	ticks(v0, v1) = collect(v0:v1)
	io = IOBuffer()
	write(io, """<rect x="$(x0)" y="$(y0)" width="$(W)" height="$(H)" fill="white" stroke="#999"/>
<text x="$(x0+W/2)" y="$(y0+15)" text-anchor="middle" font-size="13" font-family="serif" font-style="italic">$title</text>
<text x="$(x0+W/2)" y="$(y0+H-6)" text-anchor="middle" font-size="11" font-family="serif">$xlab</text>
<text x="$(x0+14)" y="$(y0+H/2)" text-anchor="middle" font-size="11" font-family="serif" transform="rotate(-90 $(x0+14) $(y0+H/2))">relerr</text>
""")
	for t in ticks(Int(xmin), Int(xmax))
		abs(xmax - xmin) < 1e-9 && break
		xx = x0 + ml + (t - xmin) / (xmax - xmin) * (W - ml - mr)
		write(io, "<line x1='$xx' y1='$(y0+mt)' x2='$xx' y2='$(y0+H-mb)' stroke='#e8e8e8'/><text x='$xx' y='$(y0+H-mb+13)' text-anchor='middle' font-size='9'>10<tspan dy='-4' font-size='7'>$t</tspan></text>\n")
	end
	for t in ticks(Int(ymin), Int(ymax))
		yy = y0 + H - mb - (t - ymin) / (ymax - ymin) * (H - mt - mb)
		write(io, "<line x1='$(x0+ml)' y1='$yy' x2='$(x0+W-mr)' y2='$yy' stroke='#e8e8e8'/><text x='$(x0+ml-5)' y='$(yy+3)' text-anchor='end' font-size='9'>10<tspan dy='-4' font-size='7'>$t</tspan></text>\n")
	end
	# 一阶/二阶参考线（锚在 vertex 首点）
	x1, x2 = xs[1], xs[end]
	y1 = ev[1]
	for (p, dash) in ((1, "2 3"), (2, "5 3"))
		y2 = y1 * (x2 / x1)^p
		(ymin <= log10(y2) <= ymax) || continue
		write(io, "<line x1='$(round(px(x1);digits=1))' y1='$(round(py(y1);digits=1))' x2='$(round(px(x2);digits=1))' y2='$(round(py(clamp(y2, 10^ymin, 10^ymax));digits=1))' stroke='#999' stroke-dasharray='$dash'/>
")
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
svg = ["<svg xmlns='http://www.w3.org/2000/svg' width='$(3W+2GAP+20)' height='$(H+52)' font-family='serif'>",
	   "<rect width='100%' height='100%' fill='white'/>",
	   "<text x='$(3W/2+GAP+10)' y='20' text-anchor='middle' font-size='15'>Toulouse IF vertex-shift vs raw: error scaling (PartialIF, single-mode ED)</text>"]
for (i, (name, xs, vs, rs, xlab)) in enumerate(panels)
	push!(svg, svg_loglog_panel(10 + (i - 1) * (W + GAP), 32, W, H, name, xs, vs, rs, xlab))
end
push!(svg, "<line x1='$(W/2-60)' y1='$(H+42)' x2='$(W/2-40)' y2='$(H+42)' stroke='#d62728' stroke-width='1.6'/><text x='$(W/2-34)' y='$(H+46)' font-size='11'>vertex</text>")
push!(svg, "<line x1='$(W/2+20)' y1='$(H+42)' x2='$(W/2+40)' y2='$(H+42)' stroke='#1f77b4' stroke-width='1.6' stroke-dasharray='7 3'/><text x='$(W/2+46)' y='$(H+46)' font-size='11'>raw (cell)</text>")
push!(svg, "<line x1='$(2W+GAP/2)' y1='$(H+42)' x2='$(2W+GAP/2+20)' y2='$(H+42)' stroke='#999' stroke-dasharray='2 3'/><text x='$(2W+GAP/2+26)' y='$(H+46)' font-size='11'>O(δ) / O(δ²) guides</text>")
push!(svg, "</svg>")
open(joinpath(OUTDIR, "vertex_shift_error_scaling.svg"), "w") do io
	write(io, join(svg, "\n"))
end
println("已保存：vertex_shift_error_scaling.svg / .csv（performance/firstorder/）")

install_vertex_accessor!()   # 恢复当前语义
install_mixed_raw_accessor!()
println("（访问器已恢复为当前语义：实时 vertex、τ 环/混合 raw）")
