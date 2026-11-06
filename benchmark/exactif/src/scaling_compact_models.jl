# =============================================================================
# src/scaling_compact_models.jl —— 严格 IF（紧凑消元路线）1/2/4 阶 Trotter 演化的
# error scaling：Toulouse 模型 + 单模 bath Anderson 模型（U≠0）
#
# 模型参数与 toulouse/、singlemodebath/ 诊断基准一致：
#   T = 0.4, β = 1.0；单模 bath = DiracDelta(ω=1.0, α=0.5)（有限模等价
#   HB=[1.0], V=[√0.5]）；Toulouse：ToulouseIM(ϵ_d=0.7), bands=1；
#   Anderson：AndersonIM(U=1.5, ϵ_d=0.7), bands=2（每带一个 bath 模）。
#
# 判据：
#   ① 机器级预检：管线 n(T) == Fock 离散 Trotter 链（同 schedule、同热态 ρ₀）
#     ——确认 U≠0 管线就是离散链的精确收缩；
#   ② error scaling：|n_N − n_exact| vs δt = T/N（n_exact = 全 H 精确对角化，
#     Toulouse 4 维 / Anderson 16 维），最小二乘阶数 p；
#   ③ 真阶判决：order 1/2 → p≈2（对角 H_imp 下严格核精确重求和 [A,B]δt²，
#     order 1 与 2 全局二阶——IFGTEMPO/-twosector 已知结论），order 4 → p≈4
#     （每减半 δt 误差 ÷16，与 order 2 的 ÷4 可区分）。
#
# 输出：results/compact_scaling_{toulouse,anderson}.{csv,png} + 终端判定表。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/scaling_compact_models.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(@__DIR__, "NativeExact.jl"))
using .NativeExact
include(joinpath(GT_ROOT, "benchmark", "exactif", "form_b", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "common.jl"))   # plot_png/PY_PLOTTER

# ---- 模型参数（诊断基准同款）----
const MODEL_T = 0.4
const BETA_M = 1.0
const EPS_D_M = 0.7
const UINT_M = 1.5
const WB, AB = 1.0, 0.5
const HB_M = ComplexF64[WB;;]
const VB_M = ComplexF64[sqrt(AB) + 0.0im]   # DiracDelta(ω=1, α=0.5) ⇔ |V|² = α

"""单带热态 e^{-βH}/tr（Toulouse）。"""
thermal2(eps_d, β) = (H = ComplexF64[0 0; 0 eps_d]; R = exp(-β .* H); R / tr(R))

"""两带热态 e^{-βH}/tr（Anderson：E = {0, ε, ε, 2ε+U}）。"""
function thermal4(eps_d, U, β)
    H = ComplexF64[0 0 0 0; 0 eps_d 0 0; 0 0 eps_d 0; 0 0 0 2eps_d + U]
    R = exp(-β .* H)
    return R / tr(R)
end

# ---------------------------------------------------------------------------
# Anderson（bands=2）16 维 Fock 参照（mode 序：imp₁=0, imp₂=1, bath₁=2, bath₂=3）
# ---------------------------------------------------------------------------

"""耦合二次型 Hq（bath ε + V 耦合，杂质能级 0）与 Him = ε(n₁+n₂) + U n₁n₂。"""
function anderson_hamiltonians()
    dim = 16
    ops = [annih_op(dim, k) for k in 0:3]
    n = [ops[k]' * ops[k] for k in 1:4]
    Hq = zeros(ComplexF64, dim, dim)
    for (imp, bath) in ((1, 3), (2, 4))
        Hq .+= HB_M[1, 1] .* n[bath]
        Hq .+= conj(VB_M[1]) .* (ops[imp]' * ops[bath])   # d†c 幅
        Hq .+= VB_M[1] .* (ops[bath]' * ops[imp])         # c†d 幅
    end
    Him = EPS_D_M .* (n[1] + n[2]) .+ UINT_M .* (n[1] * n[2])
    return Hq, Him, n
end

"""Anderson 的 bath 直积热态嵌入（bits 2 = bath₁, bits 3 = bath₂ 独立同谱）。"""
function anderson_rho_bath()
    rb = exp(-BETA_M .* HB_M[1, 1] .* ComplexF64[0 0; 0 1])
    rb = rb / tr(rb)
    r = zeros(ComplexF64, 16, 16)
    for b1 in 0:1, b2 in 0:1, b1p in 0:1, b2p in 0:1
        b = (b2 << 1 | b1)   # bath bits → 4×4 块索引
        bp = (b2p << 1 | b1p)
        for a in 0:3
            r[(b << 2 | a) + 1, (bp << 2 | a) + 1] =
                rb[b1 + 1, b1p + 1] * rb[b2 + 1, b2p + 1]
        end
    end
    return r
end

"""Anderson 参照：mode=:chain（离散链）或 :exact（连续精确）。返回 n₁(T)。"""
function anderson_ref(; mode::Symbol, order::Int=2, N::Int=1)
    Hq, Him, n = anderson_hamiltonians()
    rho0 = thermal4(EPS_D_M, UINT_M, BETA_M)
    r0 = zeros(ComplexF64, 16, 16)
    for b in 0:3, a in 0:3, a2 in 0:3
        r0[(b << 2 | a) + 1, (b << 2 | a2) + 1] = rho0[a + 1, a2 + 1]
    end
    rb = anderson_rho_bath()
    if mode === :chain
        sched = trotter_schedule(order, N, MODEL_T / N)
        U = chain_operator(Hq, Him, collect(Float64.(sched.cells)),
                           collect(Float64.(sched.links)))
    else
        U = exp(-1im * MODEL_T .* (Hq + Him))
    end
    rtot = U * (r0 * rb) * adjoint(U)
    return real(tr(rtot * n[1])) / real(tr(rtot))
end

"""Toulouse 参照（复用 formb_dense 的链机器；连续 = 4 维精确对角化）。"""
function toulouse_ref(; mode::Symbol, order::Int=2, N::Int=1)
    up = FiniteBath(HB_M, VB_M, BETA_M)
    if mode === :chain
        sched = trotter_schedule(order, N, MODEL_T / N)
        ρ = fock_chain_single(up, EPS_D_M, sched, thermal2(EPS_D_M, BETA_M))
        return real(ρ[2, 2]) / real(tr(ρ))
    end
    m = 1
    dim = 1 << (1 + m)
    ops = [annih_op(dim, k) for k in 0:m]
    hQ = bath_one_particle_hamiltonian(up)
    Hq = zeros(ComplexF64, dim, dim)
    for i in 1:(m + 1), j in 1:(m + 1)
        iszero(hQ[i, j]) || (Hq .+= hQ[i, j] .* (ops[i]' * ops[j]))
    end
    Him = EPS_D_M .* (ops[1]' * ops[1])
    Hb = zeros(ComplexF64, dim, dim)
    for i in 2:(m + 1), j in 2:(m + 1)
        iszero(HB_M[i - 1, j - 1]) || (Hb .+= HB_M[i - 1, j - 1] .* (ops[i]' * ops[j]))
    end
    Zb = tr(exp(-BETA_M .* Hb)) / 2
    rb = exp(-BETA_M .* Hb) ./ Zb
    r0 = zeros(ComplexF64, dim, dim)
    ρ0 = thermal2(EPS_D_M, BETA_M)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        r0[(b << 1 | a) + 1, (b << 1 | a2) + 1] = ρ0[a + 1, a2 + 1]
    end
    U = exp(-1im * MODEL_T .* (Hq + Him))
    rtot = U * (r0 * rb) * adjoint(U)
    ρ = zeros(ComplexF64, 2, 2)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        ρ[a + 1, a2 + 1] += rtot[(b << 1 | a) + 1, (b << 1 | a2) + 1]
    end
    return real(ρ[2, 2]) / real(tr(ρ))
end

# ---------------------------------------------------------------------------
# 紧凑管线 n(T)（折返腿 GTerm 投影读出，band 1 占据）
# ---------------------------------------------------------------------------

function compact_nt_model(model, bands::Int; order::Int, N::Int,
                          trunc=GTEMPO.truncdimcutoff(D=512, ϵ=1e-12))
    pipe = NativeExact.exactif_native_compact(HB_M, VB_M, model, N, MODEL_T / N;
                                              order=order, trunc=trunc, β=BETA_M,
                                              bands=bands, counter_c=-1)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k = lat.k
    mpsN = copy(mpsK)
    pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=1)]
    pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=1)]
    GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
    return real(GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z)
end

_slope(xs, ys) = (sum((xs .- sum(xs) / length(xs)) .* (ys .- sum(ys) / length(ys))) /
                  sum(abs2, xs .- sum(xs) / length(xs)))

"""机器级预检：管线 vs Fock 离散链（确认 U≠0 管线 = 离散链的精确收缩）。"""
function precheck()
    println("═"^78)
    println("机器级预检：紧凑管线 vs Fock 离散 Trotter 链（同 schedule、同热态 ρ₀）")
    println("═"^78)
    notrunc = NoTruncation()
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    for (tag, model, bands) in (("Toulouse", ToulouseIM(ϵ_d=EPS_D_M), 1),
                                ("Anderson", AndersonIM(U=UINT_M, ϵ_d=EPS_D_M), 2))
        for (order, N) in ((2, 1), (2, 2), (4, 1), (4, 2))
            ref = bands == 1 ?
                  toulouse_ref(mode=:chain, order=order, N=N) :
                  anderson_ref(mode=:chain, order=order, N=N)
            nt = compact_nt_model(model, bands; order=order, N=N,
                                  trunc=(N == 1 ? notrunc : trunc))
            @printf("  %-8s order %d N=%d : n = %.12f | chain %.12f | Δ = %.3e\n",
                    tag, order, N, nt, ref, abs(nt - ref))
        end
    end
end

"""标度主流程：单模型 × orders 1/2/4 × N 网格 → CSV/PNG + 判定。"""
function scaling_model(tag::String, model, bands::Int, ref_exact::Float64)
    println("═"^78)
    @printf("%s：error scaling vs 精确连续演化（n_exact = %.12f）\n", tag, ref_exact)
    println("═"^78)
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    Ns = [1, 2, 4, 8]
    csv = joinpath(GT_ROOT, "benchmark", "exactif", "results",
                   "compact_scaling_$(lowercase(tag)).csv")
    png = joinpath(GT_ROOT, "benchmark", "exactif", "results",
                   "compact_scaling_$(lowercase(tag)).png")
    names = String[]
    errs = Vector{Vector{Float64}}()
    open(csv, "w") do io
        println(io, "series,dt,relerr")
        for N in Ns
            row = Pair{String,Float64}[]
            for order in (1, 2, 4)
                nt = compact_nt_model(model, bands; order=order, N=N, trunc=trunc)
                push!(row, "exact IF (order $order)" => abs(nt - ref_exact))
            end
            if isempty(names)
                append!(names, first.(row))
                errs = [Float64[] for _ in names]
                @printf("%-6s%-10s%s\n", "N", "δt",
                        join(rpad("$(n) Δn", 16) for n in names))
            end
            @printf("%-6d%-10.6g%s\n", N, MODEL_T / N,
                    join(rpad(@sprintf("%.6e", e), 16) for e in last.(row)))
            for (j, e) in enumerate(last.(row))
                push!(errs[j], e)
                println(io, names[j], ",", MODEL_T / N, ",", e)
            end
        end
    end
    hs = log.(MODEL_T ./ Ns)
    ps = [_slope(hs, log.(errs[j])) for j in eachindex(names)]
    println("拟合阶数：", join((@sprintf("%s p = %.2f", names[j], ps[j])
                              for j in eachindex(names)), "  "))
    # 每减半 δt 的误差比（order 2 应 ≈4，order 4 应 ≈16 —— 真四阶的可区分判据）
    println("每减半 δt 误差比：")
    for (j, name) in enumerate(names)
        ratios = [errs[j][i - 1] / errs[j][i] for i in 2:length(Ns)]
        println("  ", name, " : ", join((@sprintf("%.1f", r) for r in ratios), ", "))
    end
    ok = abs(ps[1] - 2) < 0.2 && abs(ps[2] - 2) < 0.2 && abs(ps[3] - 4) < 0.4
    println("真阶判定（1/2 → p≈2，4 → p≈4）：", ok ? "PASS" : "FAIL")
    title = "strict-IF (exact kernel) Trotter error scaling — $tag\n" *
            "T=$(MODEL_T), β=$(BETA_M), single-mode bath (ω=$WB, α=$AB)" *
            (tag == "Anderson" ? ", U=$(UINT_M), ϵ_d=$(EPS_D_M), bands=2" :
             ", ϵ_d=$(EPS_D_M)")
    try
        plot_png(csv, png, title, "δt = T / N", "|n(T) − n_exact| (exact-IF compact route)")
        println("图 -> ", png)
    catch e
        println("画图失败（CSV 已写）：", sprint(showerror, e))
    end
    return ps
end

function main()
    precheck()
    ps_t = scaling_model("Toulouse", ToulouseIM(ϵ_d=EPS_D_M), 1,
                         toulouse_ref(mode=:exact))
    ps_a = scaling_model("Anderson", AndersonIM(U=UINT_M, ϵ_d=EPS_D_M), 2,
                         anderson_ref(mode=:exact))
    println("═"^78)
    @printf("总结：Toulouse p = %.2f/%.2f/%.2f | Anderson p = %.2f/%.2f/%.2f（orders 1/2/4）\n",
            ps_t..., ps_a...)
end

main()
