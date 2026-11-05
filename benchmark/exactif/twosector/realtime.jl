# =============================================================================
# exactif/twosector/realtime.jl —— 两扇区严格核基准（与 IFGTEMPO 对齐）
#
# 单模/有限 bath 上，严格核 exact_Δt 的 1/2/4 阶 Trotter 演化的权威标度：
# 采用 IFGTEMPO 的两扇区构造（IF 变量与 IMP 轨迹变量分离，接口测度
# M = exp[−ā_IF x_K − x̄_K a_IF] 交叉耦合两侧——见 IFGTEMPO/docs/
# IF_DERIVATION_AND_PACKAGE_LOGIC.md 第 3、9 节），对精确连续 Fock 演化
# 验证密度矩阵误差的收敛阶数：
#
#   order 1（全步长放置）      → 全局 O(δt²)（精确顶点幅重求和 [A,B]δt² 首项）
#   order 2（Strang 半全半）   → 全局 O(δt²)
#   order 4（Suzuki S2²S2）    → 全局 O(δt⁴)
#
# 该脚本依赖外部 IFGTEMPO 包（GTEMPO 的姊妹仓库，path 引入；GTEMPO 的
# GrassmannMPS 机制是其全部底层实现）。单扇区 GTEMPO 管线（../common.jl）
# 与此构造的差异见 benchmark/exactif/README.md 的"与 GTEMPO 的差异"一节。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/twosector/realtime.jl
# =============================================================================
using LinearAlgebra, Printf
using GTEMPO

# ---- 引入 IFGTEMPO（姊妹仓库，与 GTEMPO 同级）----
const IFT_ROOT = joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO")
isdir(IFT_ROOT) || error("IFGTEMPO 仓库未找到：$IFT_ROOT（两扇区基准需要它）")
push!(LOAD_PATH, IFT_ROOT)
println("引入 IFGTEMPO：", IFT_ROOT)
using IFGTEMPO

const MODEL_T = 0.34      # 总演化时长
const D_TRUNC = typemax(Int)
const EPS_TRUNC = 1.0e-12

"""两自旋 + M 模式有限 bath 的构造（每自旋独立 bath）。"""
function two_spin_setup(M::Int)
    M == 1 && return (FiniteBath([0.4;;], [0.23 + 0.11im], 1.3),
                      FiniteBath([0.35;;], [0.19 - 0.07im], 1.3))
    M == 2 && return (FiniteBath([0.4 0.1; 0.1 0.3], [0.23 + 0.11im, 0.15 - 0.05im], 1.3),
                      FiniteBath([0.38 0.08; 0.08 0.28], [0.19 - 0.07im, 0.12 + 0.04im], 1.3))
    error("M = $M 未预置")
end

"""精确连续 Fock 演化的约化密度矩阵 ρ_imp(T)。"""
function exact_density(up, dn, imp, T::Real)
    mu, md = bath_modes(up), bath_modes(dn)
    btot = mu + md; io = btot; total = btot + 2; dim = 1 << total
    ops = [annihilation_operator(total, mode) for mode in 0:(total - 1)]
    Hq = zeros(ComplexF64, dim, dim)
    for (mdl, off, impmode) in ((up, 0, io), (dn, mu, io + 1))
        for i in 1:bath_modes(mdl), j in 1:bath_modes(mdl)
            Hq .+= mdl.bath_hamiltonian[i, j] .* (adjoint(ops[off + i]) * ops[off + j])
        end
        for i in 1:bath_modes(mdl)
            Hq .+= conj(mdl.hybridization[i]) .* (adjoint(ops[impmode + 1]) * ops[off + i])
            Hq .+= mdl.hybridization[i] .* (adjoint(ops[off + i]) * ops[impmode + 1])
        end
    end
    Himp = imp.epsilon_up .* number_operator(total, io) +
           imp.epsilon_down .* number_operator(total, io + 1) +
           imp.interaction_u .* (number_operator(total, io) * number_operator(total, io + 1))
    U = exp((-1im * T) .* (Hq + Himp))
    rho_b = kron(thermal_bath_state(dn), thermal_bath_state(up))
    rho = U * kron(imp.initial_density, rho_b) * adjoint(U)
    reduced = zeros(ComplexF64, 4, 4); bath_dim = 1 << btot
    for i in 0:3, j in 0:3, bs in 0:(bath_dim - 1)
        reduced[i + 1, j + 1] += rho[bs + (i << btot) + 1, bs + (j << btot) + 1]
    end
    return reduced
end

function path_density_err(order::Int, N::Int, up, dn, imp, T::Real)
    path = build_impurity_path(up, dn, N, T / N; order = order,
                               impurity_model = imp, cutoff = EPS_TRUNC)
    rho = gtempo_density_matrix_direct(path.factors...)
    return maximum(abs.(rho .- exact_density(up, dn, imp, T)))
end

function main(; M::Int=2, Ns=(2, 4, 8), T::Real=MODEL_T,
               eps_up::Real=0.35, eps_dn::Real=-0.22, Uint::Real=1.0, label::String="")
    up, dn = two_spin_setup(M)
    rho0 = ComplexF64[0 0 0 0; 0 0.3 0 0; 0 0 0.2 0; 0 0 0 0.5]
    imp = TwoSpinImpurityModel(eps_up, eps_dn, Uint, rho0)

    println("=" ^ 100)
    @printf("两扇区严格核（IFGTEMPO 对齐）Trotter 标度 —— 两自旋 (ε=%.2f,%.2f, U=%.2f)，M=%d 模式 bath，T=%.2f\n",
            eps_up, eps_dn, Uint, M, T)
    println("=" ^ 100)
    csv = joinpath(@__DIR__, "realtime.csv")
    results = Tuple{String,String,Float64,Float64}[]   # (model, order, dt, err)
    open(csv, "w") do io
        println(io, "series,order,dt,err")  # series = "<model> | order k"
        for order in (1, 2, 4)
            errs = Tuple{Float64,Float64}[]
            for N in Ns
                t0 = time()
                e = path_density_err(order, N, up, dn, imp, T)
                push!(errs, (T / N, e))
                println(io, replace(label, "," => ";"), " | order $(order)", ",", N, ",", T / N, ",", e)
                push!(results, (label, "order $(order)", T / N, e))
                @printf("  order %d N=%2d: err = %.6e  (%.1fs)\n", order, N, e, time() - t0)
            end
            p1 = log(errs[1][2] / errs[2][2]) / log(2)
            p2 = length(Ns) > 2 ? log(errs[2][2] / errs[3][2]) / log(2) : NaN
            @printf("  order %d 汇总: err(δt)=%.3e  err(δt/2)=%.3e  err(δt/4)=%.3e   p = %.2f, %.2f\n\n",
                    order, errs[1][2], errs[2][2], errs[3][2], p1, p2)
        end
    end
    println("数据 -> ", csv)
    png = joinpath(@__DIR__, "realtime_" * replace(label, "," => ";", " " => "") * ".png")
    plot_orders(results, png, T)
    println("图   -> ", png)
    return nothing
end

# log-log 标度图：每个 (model, order) 一条曲线，标注最小二乘拟合阶；
# O(δt²)/O(δt⁴) 参考线锚在各自首条曲线首点。
const PY_PLOTTER = raw"""
import sys, csv, math
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path, png_path = sys.argv[1:3]
series = defaultdict(list)
models = []
with open(csv_path) as f:
    for row in csv.DictReader(f):
        series[row["series"]].append((float(row["dt"]), float(row["err"])))
        mdl = row["series"].split(" | ")[0]
        if mdl not in models:
            models.append(mdl)

colors = {"order 1": "#7f7f7f", "order 2": "#1f77b4", "order 4": "#d62728"}
fig, axes = plt.subplots(1, len(models), figsize=(6.2 * len(models), 5.0), sharey=True)
if len(models) == 1:
    axes = [axes]
for ax, mdl in zip(axes, models):
    for od in ("order 1", "order 2", "order 4"):
        key = mdl + " | " + od
        if key not in series:
            continue
        pts = sorted(series[key])
        dt = [q[0] for q in pts]
        er = [q[1] for q in pts]
        ax.loglog(dt, er, "o-", color=colors[od], label=od + " placement")
        xx = [math.log(v) for v in dt]
        yy = [math.log(v) for v in er]
        n = len(xx)
        xm, ym = sum(xx) / n, sum(yy) / n
        p = sum((a - xm) * (b - ym) for a, b in zip(xx, yy)) / sum((a - xm) ** 2 for a in xx)
        ax.annotate("p = {0:.2f}".format(p), xy=(dt[-1], er[-1]), xytext=(5, 4),
                    textcoords="offset points", fontsize=10, color=colors[od])
    # reference slopes anchored at the order-2 curve head
    key2 = mdl + " | order 2"
    pts = sorted(series[key2])
    (d0, e0) = pts[0]
    d2 = pts[-1][0]
    for pexp, ls in ((2, ":"), (4, "-.")):
        ax.loglog([d0, d2 * 1.15], [e0, e0 * (d2 * 1.15 / d0) ** pexp], ls=ls, color="gray", lw=1)
        ax.annotate(r"$O(\delta t^{0})$".format(pexp),
                    xy=(d2 * 1.15, e0 * (d2 * 1.15 / d0) ** pexp),
                    xytext=(-6, 8), textcoords="offset points", fontsize=10, color="gray")
    ax.set_title(mdl, fontsize=11)
    ax.set_xlabel(r"$\delta t = T / N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
axes[0].set_ylabel(r"max $|\rho_{imp} - \rho_{exact}|$")
fig.suptitle("Two-sector exact-IF kernel (aligned with IFGTEMPO): Trotter-order scaling\n"
             "vs the exact continuous Fock evolution", fontsize=12)
fig.tight_layout()
fig.savefig(png_path, dpi=150)
print("png ->", png_path)
"""

function plot_orders(results, png::String, T::Real)
    py = tempname()
    open(py, "w") do io
        write(io, PY_PLOTTER)
    end
    csv = tempname()
    open(csv, "w") do io
        println(io, "series,dt,err")
        for (mdl, od, dt, err) in results
            println(io, replace(mdl, "," => ";"), " | ", od, ",", dt, ",", err)
        end
    end
    try
        run(pipeline(`python3 $py $csv $png`; stdout=stdout, stderr=stderr))
    finally
        rm(py; force=true)
        rm(csv; force=true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(label="two-spin (0.35,-0.22,U=1.0)")
    main(label="Anderson-like (0.7,0.7,U=1.5)", eps_up=0.7, eps_dn=0.7, Uint=1.5)
end
