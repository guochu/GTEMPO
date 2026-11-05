# =============================================================================
# src/scaling_formb.jl —— FormB（消元式单扇区 GMPS）Trotter 阶数标度
#
# 对精确连续 Fock 演化的密度矩阵误差 vs δt（log-log）：
#   order 1/2（精确顶点幅已重求和 [A,B]δt² 首项）→ 全局二阶
#   order 4（Suzuki）→ 四阶
# 叠加 IFGTEMPO 两扇区 GMPS（build_impurity_path，同一离散链）的参照点
# ——两条曲线应重合（FormB 与 IFGTEMPO 对齐的直接可视化）。
#
# 内存纪律：一律带截断 truncdimcutoff(D=512, ϵ=1e-12)（Δ_eff 非局域 →
# bond 按 4^N 增长，必须封顶）；单自旋/两自旋分两个 julia 进程跑
# （ARGS = single | twospin | plot）。
#
# 输出：benchmark/exactif/results/formb_scaling_<model>.csv / .png
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/scaling_formb.jl single
#       julia --project=<GTEMPO 根> benchmark/exactif/src/scaling_formb.jl twospin
#       julia --project=<GTEMPO 根> benchmark/exactif/src/scaling_formb.jl plot
# =============================================================================
using LinearAlgebra, Printf
using GTEMPO

const IFT_ROOT = joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO")
isdir(IFT_ROOT) || error("IFGTEMPO 仓库未找到：$IFT_ROOT")
push!(LOAD_PATH, IFT_ROOT)
using IFGTEMPO

include(joinpath(@__DIR__, "..", "form_b", "formb_dense.jl"))   # Fock 参照件
include(joinpath(@__DIR__, "FormB.jl"))
using .FormB

const T_TOT = 0.4
const BETA = 1.3
const EPS_D = 0.7
const RESULTS = joinpath(@__DIR__, "..", "results")
isdir(RESULTS) || mkdir(RESULTS)

const TRUNC = GTEMPO.truncdimcutoff(D = 512, ϵ = 1e-12, add_back = 0)

# ---------------------------------------------------------------------------
# 精确连续演化参照（ρ_exact = Tr_bath[e^{-iT H} ρ₀⊗ρ_B e^{+iT H}]）
# ---------------------------------------------------------------------------

"""单自旋连续演化（mode 序 imp=0, bath≥1）。"""
function exact_single(up::FiniteBath, eps_d::Real, T::Real, rho0::AbstractMatrix)
    m = bath_modes(up)
    nt = 1 + m
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt - 1)]
    hQ = bath_one_particle_hamiltonian(up)
    Hq = zeros(ComplexF64, dim, dim)
    for i in 1:(m + 1), j in 1:(m + 1)
        iszero(hQ[i, j]) && continue
        Hq .+= hQ[i, j] .* (ops[i]' * ops[j])
    end
    Htot = Hq + eps_d .* (ops[1]' * ops[1])
    Hb = zeros(ComplexF64, dim, dim)
    for i in 2:(m + 1), j in 2:(m + 1)
        iszero(up.bath_hamiltonian[i - 1, j - 1]) && continue
        Hb .+= up.bath_hamiltonian[i - 1, j - 1] .* (ops[i]' * ops[j])
    end
    Zb = tr(exp(-up.beta .* Hb)) / 2
    rb = exp(-up.beta .* Hb) ./ Zb
    r0 = zeros(ComplexF64, dim, dim)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        r0[(b << 1 | a) + 1, (b << 1 | a2) + 1] = rho0[a + 1, a2 + 1]
    end
    U = exp(-1im * T .* Htot)
    rtot = U * (r0 * rb) * adjoint(U)
    rho = zeros(ComplexF64, 2, 2)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        rho[a + 1, a2 + 1] += rtot[(b << 1 | a) + 1, (b << 1 | a2) + 1]
    end
    return rho
end

"""两自旋连续演化（mode 序 imp_up=0, imp_dn=1, up-bath, dn-bath）。"""
function exact_two_spin(up::FiniteBath, dn::FiniteBath,
                        imp::TwoSpinImpurityModel, T::Real)
    mu, md = bath_modes(up), bath_modes(dn)
    nt = 2 + mu + md
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt - 1)]
    Hq = zeros(ComplexF64, dim, dim)
    Hb = zeros(ComplexF64, dim, dim)
    for (model, off, impmode) in ((up, 2, 0), (dn, 2 + mu, 1))
        for i in 1:bath_modes(model), j in 1:bath_modes(model)
            Hq .+= model.bath_hamiltonian[i, j] .* (ops[off + i]' * ops[off + j])
            Hb .+= model.bath_hamiltonian[i, j] .* (ops[off + i]' * ops[off + j])
        end
        for i in 1:bath_modes(model)
            Hq .+= conj(model.hybridization[i]) .* (ops[impmode + 1]' * ops[off + i])
            Hq .+= model.hybridization[i] .* (ops[off + i]' * ops[impmode + 1])
        end
    end
    Him = imp.epsilon_up .* (ops[1]' * ops[1]) +
          imp.epsilon_down .* (ops[2]' * ops[2]) +
          imp.interaction_u .* (ops[1]' * ops[1]) * (ops[2]' * ops[2])
    Zb = tr(exp(-up.beta .* Hb)) / 4
    rb = exp(-up.beta .* Hb) ./ Zb
    r0 = zeros(ComplexF64, dim, dim)
    for b in 0:((1 << (nt - 2)) - 1), a in 0:3, a2 in 0:3
        r0[(b << 2 | a) + 1, (b << 2 | a2) + 1] = imp.initial_density[a + 1, a2 + 1]
    end
    U = exp(-1im * T .* (Hq + Him))
    rtot = U * (r0 * rb) * adjoint(U)
    rho = zeros(ComplexF64, 4, 4)
    for b in 0:((1 << (nt - 2)) - 1), a in 0:3, a2 in 0:3
        rho[a + 1, a2 + 1] += rtot[(b << 2 | a) + 1, (b << 2 | a2) + 1]
    end
    return rho
end

# ---------------------------------------------------------------------------
# 标度扫描：每 (order, N) 记录 FormB 与 IFGTEMPO 两扇区的误差
# ---------------------------------------------------------------------------

"""追加一行到 CSV（不存在则写表头）。"""
function append_row(csv::String, model::String, order::Int, N::Int, dt::Real,
                    err_formb::Real, err_ifg::Real)
    header = !isfile(csv)
    open(csv, "a") do io
        header && println(io, "model,order,N,dt,err_formb,err_ifgtempo")
        println(io, model, ",", order, ",", N, ",", dt, ",", err_formb, ",", err_ifg)
    end
end

function run_scaling_single()
    up = FiniteBath([0.4;;], [0.23 + 0.11im], BETA)
    rho0 = ComplexF64[1 0; 0 0]
    model = "single-spin (eps=0.7; M=1)"
    rho_ex = exact_single(up, EPS_D, T_TOT, rho0)
    csv = joinpath(RESULTS, "formb_scaling_single.csv")
    isfile(csv) && rm(csv)
    println("单自旋标度（T=$(T_TOT)，D≤512 截断）：参照 tr ρ = ", real(tr(rho_ex)))
    for order in (1, 2, 4)
        for N in (2, 4, 8, 16)
            dt = T_TOT / N
            # FormB（消元式单扇区，GTEMPO GMPS）
            path = FormB.formb_path((up,), (eps = EPS_D, rho0 = rho0), N, dt;
                                    order = order, counter_c = -1, trunc = TRUNC)
            X = path.factors
            rho_fb = FormB.formb_density_matrix_single(X...)
            err_fb = maximum(abs.(rho_fb - rho_ex))
            # IFGTEMPO 两扇区参照（同一离散链：下自旋完全解耦 V=0, ε=0，
            # ρ₀ = |0↑⟩⟨0↑| ⊗ |0↓⟩⟨0↓|；上自旋约化 = ρ₄[1:2, 1:2]）
            dn = FiniteBath([0.0;;], [0.0 + 0im], BETA)
            rho0ker = ComplexF64[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
            imp2 = TwoSpinImpurityModel(EPS_D, 0.0, 0.0;
                                        initial_density = rho0ker)
            path2 = build_impurity_path(up, dn, N, dt; order = order,
                                        impurity_model = imp2, cutoff = 1e-12)
            rho_if = gtempo_density_matrix_direct(path2.factors...)
            err_if = maximum(abs.(rho_if[1:2, 1:2] - rho_ex))
            append_row(csv, model, order, N, dt, err_fb, err_if)
            @printf("  order %d N=%2d : dt=%.5f  err_formb=%.3e  err_ifgtempo=%.3e\n",
                    order, N, dt, err_fb, err_if)
            GC.gc()
        end
    end
    println("CSV -> ", csv)
end

function run_scaling_twospin()
    # T = 0.34：与 twosector 发表基准/IFGTEMPO test_trotter_orders 同参数
    #（T=0.4 时 order 4 进入固定总时间平台，比值 ~10）
    T2 = 0.34
    up = FiniteBath([0.4 0.1; 0.1 0.3], [0.23 + 0.11im, 0.15 - 0.05im], BETA)
    dn = FiniteBath([0.38 0.08; 0.08 0.28], [0.19 - 0.07im, 0.12 + 0.04im], BETA)
    # ρ₀ = diag(0, 0.3, 0.2, 0.5)：与 twosector 发表基准完全同初态
    #（有占据的对角初态；order 4 平台行为随初态不同）
    rho0_2 = ComplexF64[0 0 0 0; 0 0.3 0 0; 0 0 0.2 0; 0 0 0 0.5]
    imp = TwoSpinImpurityModel(0.35, -0.22, 1.0; initial_density = rho0_2)
    model = "two-spin (0.35;-0.22;U=1.0;T=0.34)"
    rho_ex = exact_two_spin(up, dn, imp, T2)
    csv = joinpath(RESULTS, "formb_scaling_twospin.csv")
    isfile(csv) && rm(csv)
    println("两自旋标度（T=$(T_TOT)，D≤512 截断）：参照 tr ρ = ", real(tr(rho_ex)))
    for order in (1, 2, 4)
        for N in (2, 4, 8, 16)
            dt = T2 / N
            # FormB
            path = FormB.formb_path((up, dn), imp, N, dt;
                                    order = order, counter_c = -1, trunc = TRUNC)
            X = path.factors
            rho_fb = FormB.formb_density_matrix(X...)
            err_fb = maximum(abs.(rho_fb - rho_ex))
            # IFGTEMPO 两扇区
            path2 = build_impurity_path(up, dn, N, dt; order = order,
                                        impurity_model = imp, cutoff = 1e-12)
            rho_if = gtempo_density_matrix_direct(path2.factors...)
            err_if = maximum(abs.(rho_if - rho_ex))
            append_row(csv, model, order, N, dt, err_fb, err_if)
            @printf("  order %d N=%2d : dt=%.5f  err_formb=%.3e  err_ifgtempo=%.3e\n",
                    order, N, dt, err_fb, err_if)
            GC.gc()
        end
    end
    println("CSV -> ", csv)
end

# ---------------------------------------------------------------------------
# 出图（CSV → matplotlib，同 twosector 的样式；FormB 实线、IFGTEMPO 空心点）
# ---------------------------------------------------------------------------

const PY_PLOTTER = raw"""
import sys, csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path, png_path, title = sys.argv[1:4]
series = defaultdict(list)
with open(csv_path) as f:
    for row in csv.DictReader(f):
        series[row["order"]].append((float(row["dt"]), float(row["err_formb"]),
                                     float(row["err_ifgtempo"])))

colors = {"1": "#7f7f7f", "2": "#1f77b4", "4": "#d62728"}
fig, ax = plt.subplots(figsize=(7.0, 5.2))
for od in ("1", "2", "4"):
    if od not in series:
        continue
    pts = sorted(series[od])
    dt   = [q[0] for q in pts]
    ef   = [q[1] for q in pts]
    ei   = [q[2] for q in pts]
    lbl = "order " + od
    ax.loglog(dt, ef, "o-", color=colors[od], label=lbl + " (FormB)")
    ax.loglog(dt, ei, "o", mfc="none", color=colors[od], label=lbl + " (IFGTEMPO)")
    xx = [math.log(v) for v in dt]
    yy = [math.log(v) for v in ef]
    n = len(xx)
    xm, ym = sum(xx) / n, sum(yy) / n
    p = sum((a - xm) * (b - ym) for a, b in zip(xx, yy)) / sum((a - xm) ** 2 for a in xx)
    ax.annotate("p = {0:.2f}".format(p), xy=(dt[-1], ef[-1]), xytext=(5, 4),
                textcoords="offset points", fontsize=10, color=colors[od])
    # 参考斜率锚在曲线首点
    d0, e0 = dt[0], ef[0]
    d1 = dt[-1]
    for pexp, ls in ((2, ":"), (4, "-.")):
        if str(pexp) == od or (pexp == 2 and od in ("1", "2")):
            ax.loglog([d0, d1 * 1.2], [e0, e0 * (d1 * 1.2 / d0) ** pexp],
                      ls=ls, color="gray", lw=1)
            ax.annotate(r"$O(\delta t^{0})$".format(pexp),
                        xy=(d1 * 1.2, e0 * (d1 * 1.2 / d0) ** pexp),
                        xytext=(-6, 8), textcoords="offset points",
                        fontsize=10, color="gray")
ax.set_xlabel(r"$\delta t = T / N$")
ax.set_ylabel(r"max $|\rho_{imp} - \rho_{exact}|$")
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="lower left", fontsize=9)
ax.set_title(title + "\nvs the exact continuous Fock evolution", fontsize=11)
fig.tight_layout()
fig.savefig(png_path, dpi=150)
print("png ->", png_path)
"""

function make_plot(csv::String, png::String, title::String)
    code = tempname() * ".py"
    open(code, "w") do io
        print(io, "import math\n", PY_PLOTTER)
    end
    run(`python3 $code $csv $png $title`)
    rm(code)
end

stage = isempty(ARGS) ? "all" : ARGS[1]
if stage in ("all", "single")
    run_scaling_single()
    make_plot(joinpath(RESULTS, "formb_scaling_single.csv"),
              joinpath(RESULTS, "formb_scaling_single.png"),
              "Form B eliminated single-sector (GTEMPO GMPS), single spin, order 1/2/4")
end
if stage in ("all", "twospin")
    run_scaling_twospin()
    make_plot(joinpath(RESULTS, "formb_scaling_twospin.csv"),
              joinpath(RESULTS, "formb_scaling_twospin.png"),
              "Form B eliminated single-sector (GTEMPO GMPS), two spin, order 1/2/4")
end
