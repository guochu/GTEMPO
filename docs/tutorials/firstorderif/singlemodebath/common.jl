# =============================================================================
# firstorderif/singlemodebath/common.jl —— 一阶影响泛函误差标度基准公共部分
#
# 模型：单模 bath 的 Anderson 杂质模型（单轨道双带，每带一个单模费米 bath）
#   H = ϵ_d (n₁ + n₂) + U n₁n₂
#       + Σ_m [ω b†_m b_m + √α (a†_m b_m + b†_m a_m)]。
# 谱密度取 DiracDelta(ω, α)：关联函数单指数、IF 键维数小，截断可取得足够紧
#（D=200, ϵ=1e-12），使测到的误差是离散化误差而非键截断误差。
#
# 各时间轴脚本（realtime.jl / imagtime.jl / …）include 本文件，只保留各自的
# GTEMPO 管线与参考解；扫描、误差表、CSV 与 log-log 标度图由 scaling_main
# 统一驱动。参考解 = 全模型（16 维 Fock 空间）的精确对角化（ED，见
# test/util.jl 的 singlemode_ed / greater_lesser_ed / gτ_ed）。
# =============================================================================
using LinearAlgebra, Printf
using GTEMPO
using Z2Tensors: TruncationScheme
using ImpurityModelBase
# test/util.jl 的 ordering 表引用了两个未导出的 ordering 结构体
using GTEMPO: A1B1B̄1Ā1, A1B1ā1b̄1Ā1B̄1a1b1

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "test", "util.jl"))
include(joinpath(GT_ROOT, "test", "normalbath", "util.jl"))

# ------------------------------------------------------------------ 模型 --
const β = 1.0            # 逆温（虚轴长度 / 热态）
const U = 1.0            # Hubbard 相互作用
const ϵ_d = 0.7          # 杂质能级
const ωb, αb = 1.0, 0.5  # 单模 bath 频率与耦合（每带一个）

# ------------------------------------------------------ 误差、阶数与输出 --
_m(v) = sum(v) / length(v)
# 误差取自然对数、横轴取 -log2(N)（∝ log δ）：斜率即「误差对 δ 的阶」
slope(hs, es) = (sum((hs .- _m(hs)) .* log.(es)) / sum(abs2, hs .- _m(hs))) / log(2)
orders(errs) = [log2(errs[i-1] / errs[i]) for i in 2:length(errs)]

# 与具体基准无关的 matplotlib 画图脚本：argv = csv, png, title, xlabel, ylabel。
# csv 列：series,dt,relerr；每个 series 画一条 log-log 线并标注最小二乘拟合阶，
# 另画锚在首条曲线首点的 O(δ)/O(δ²) 参考线。
const PY_PLOTTER = raw"""
import sys, csv, math
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path, png_path, title, xlabel, ylabel = sys.argv[1:6]
series = defaultdict(list)
with open(csv_path) as f:
    for row in csv.DictReader(f):
        series[row["series"]].append((float(row["dt"]), float(row["relerr"])))

fig, ax = plt.subplots(figsize=(7.0, 5.0))
colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
dts = None
for i, ((name, pts), c) in enumerate(zip(sorted(series.items()), colors)):
    pts.sort()
    dt = [p[0] for p in pts]
    er = [p[1] for p in pts]
    dts = dt
    ax.loglog(dt, er, "o-", color=c, label=name)
    x = [math.log(d) for d in dt]
    y = [math.log(e) for e in er]
    n = len(x)
    xm, ym = sum(x) / n, sum(y) / n
    p = sum((xi - xm) * (yi - ym) for xi, yi in zip(x, y)) / \
        sum((xi - xm) ** 2 for xi in x)
    # 相邻 series 上下错开标注，避免曲线重合时标注互相遮挡
    ax.annotate(f"order = {p:.2f}", xy=(dt[-1], er[-1]),
                xytext=(4, -14 if i % 2 == 0 else 10),
                textcoords="offset points", fontsize=9, color=c)

# O(dt) / O(dt^2) 参考线：锚在第一条曲线的首点
_, pts0 = sorted(series.items())[0]
d1, e1 = min(pts0)
d2 = max(dts)
for p, ls in ((1, "--"), (2, ":")):
    ax.loglog([d1, d2], [e1, e1 * (d2 / d1) ** p], ls, color="gray", lw=1,
              label=rf"$O(\delta t^{p})$")

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(title, fontsize=11)
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="best", fontsize=9)
fig.tight_layout()
fig.savefig(png_path, dpi=150)
print("图 ->", png_path)
"""

function plot_png(csv, png, title, xlabel, ylabel)
    py = tempname()
    open(py, "w") do io
        write(io, PY_PLOTTER)
    end
    try
        run(`python3 $py $csv $png $title $xlabel $ylabel`)
    finally
        rm(py; force=true)
    end
end

"""
    scaling_main(runner; Ns, dtof, nlabel, csvfile, pngfile, title, xlabel, ylabel, D, ϵ)

通用标度驱动。`runner(N; trunc)` 跑一条管线，返回 `Pair("序列名", 相对误差)`
的有序向量（如 `[("greater gt", e1), ("lesser lt", e2)]`）。扫描 Ns：打印误差表
（含相邻误差比的阶数）、写 CSV（series,dt,relerr）、拟合 log-log 斜率并画标度图。
"""
function scaling_main(runner; Ns, dtof, nlabel="N", csvfile, pngfile, title,
                      xlabel, ylabel, D=200, ϵ=1.0e-12)
    trunc = truncdimcutoff(D=D, ϵ=ϵ)
    hs = -log2.(collect(Ns))
    names = String[]
    errs = Vector{Vector{Float64}}()   # errs[j][i]：第 j 条序列在 Ns[i] 的误差
    open(csvfile, "w") do io
        println(io, "series,dt,relerr")
        for (i, N) in enumerate(Ns)
            t0 = time()
            res = runner(N; trunc=trunc)
            if isempty(names)
                append!(names, first.(res))
                errs = [Float64[] for _ in names]
                @printf("%-6s%-12s%s%s\n", nlabel, "δ",
                        join(rpad("$(names[j]) relerr", 22) for j in eachindex(names)),
                        join((j > 1 ? "| " : "") * rpad(names[j], 6)
                             for j in eachindex(names)))
            end
            es = last.(res)
            for (j, e) in enumerate(es)
                push!(errs[j], e)
                println(io, names[j], ",", dtof(N), ",", e)
            end
            ostr = join((j > 1 ? "| " : "") *
                        (i == 1 ? rpad("  --", 6) :
                         rpad(@sprintf("%4.2f", log2(errs[j][i-1] / es[j])), 6))
                        for j in eachindex(es))
            @printf("%-6d%-12.6g%s%s  (%.1fs)\n", N, dtof(N),
                    join(rpad(@sprintf("%.6e", e), 22) for e in es), ostr, time() - t0)
        end
    end
    println("拟合收敛阶：",
            join((rpad("$(names[j]) p = $(@sprintf("%.2f", slope(hs, errs[j])))", 24)
                  for j in eachindex(names))), "（一阶 scaling 即 p ≈ 1）")
    plot_png(csvfile, pngfile, title, xlabel, ylabel)
    println("数据 -> ", csvfile)
    return errs
end
