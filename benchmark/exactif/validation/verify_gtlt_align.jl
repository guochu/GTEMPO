# =============================================================================
# validation/verify_gtlt_align.jl —— 严格核 gt/lt Green 函数系列：
# 原生紧凑管线 vs IFGTEMPO 的机器级对齐（①单带 + ⑤双带 U=1.5）
#
# 约定（IFGTEMPO gtempo_greens_function_direct，已测试覆盖）：
#   c(t_f) @ imp(N+1,−,eta)；greater[k+1] = ⟨c(t_f)c†(k·δt)⟩（fwd_bar@k+1）；
#   lesser[k+1] = ⟨c†(k·δt)c(t_f)⟩（bwd_bar@k）；Z 归一；裸关联。
# 原生镜像 = NativeExact.gtlt_compact（节点 = 层 + 1）。
#
# 判据：
#   ① 单带 Toulouse（order 2 N=3，order 4 N=1）：native vs IFGTEMPO 解耦双自旋
#     flavor-1，机器级；
#   ⑤ 双带 Anderson（U=1.5，热态 4×4 初态，order 2 N=3）：native bands=2 vs
#     IFGTEMPO 直接同模型 flavor=1/2，机器级；
#   附：native gt/lt 对 ED 连续演化在有效时刻的 Trotter 收敛（dt 减半误差下降）。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/verify_gtlt_align.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(@__DIR__, "..", "src", "NativeExact.jl"))
using .NativeExact
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
using GTEMPO: A1B1B̄1Ā1, A1B1ā1b̄1Ā1B̄1a1b1
include(joinpath(GT_ROOT, "test", "util.jl"))
include(joinpath(GT_ROOT, "test", "normalbath", "util.jl"))

push!(LOAD_PATH, joinpath(dirname(dirname(dirname(pathof(GTEMPO)))), "IFGTEMPO"))
using IFGTEMPO

const T = 0.4
const BETA = 1.0
const HB1 = 1.0
const VB1 = sqrt(0.5)
const EPS_D = 0.7
const UINT = 1.5
const HB_M = ComplexF64[HB1;;]
const VB_M = ComplexF64[VB1 + 0.0im]

thermal4(ε, U, β) = (H = ComplexF64[0 0 0 0; 0 ε 0 0; 0 0 ε 0; 0 0 0 2ε + U];
                     R = exp(-β .* H); R / tr(R))
thermal2(ε, β) = (H = ComplexF64[0 0; 0 ε]; R = exp(-β .* H); R / tr(R))

"""原生紧凑管线的 gt/lt（order, N, bands, U 或 0）。"""
function native_gtlt(; order::Int, N::Int, bands::Int, U::Real, eps_d::Real=EPS_D,
                     β::Real=BETA, ρ0::Union{Nothing,Matrix{ComplexF64}}=nothing)
    model = U == 0 ? ToulouseIM(ϵ_d=eps_d) : AndersonIM(U=U, ϵ_d=eps_d)
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    pipe = NativeExact.exactif_native_compact(HB_M, VB_M, model, N, T / N;
                                              order=order, trunc=trunc, β=β,
                                              bands=bands, counter_c=-1)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k_node = lat.k
    csite = GrassmannContourIndex(k_node, conj=false, branch=:-, band=1)
    gt = zeros(ComplexF64, length(pipe.cells) + 1)
    lt = zeros(ComplexF64, length(pipe.cells) + 1)
    for k in 0:length(pipe.cells)
        g1 = GrassmannContourIndex(k + 2, conj=true, branch=:+, band=1)
        mpsN = copy(mpsK)
        GTEMPO.apply!(GTEMPO.GTerm(lat[g1], lat[csite], coeff=1.0), mpsN)
        gt[k + 1] = GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z
        l1 = GrassmannContourIndex(k + 1, conj=true, branch=:-, band=1)
        mpsN = copy(mpsK)
        GTEMPO.apply!(GTEMPO.GTerm(lat[l1], lat[csite], coeff=1.0), mpsN)
        lt[k + 1] = GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z
    end
    return gt, lt
end

"""IFGTEMPO 两自旋路线的 gt/lt（U=0 解耦或 U=1.5 同模型；flavor=1 即 band-1）。"""
function ifgtempo_gtlt(; order::Int, N::Int, U::Real, eps_d::Real=EPS_D,
                       β::Real=BETA, ρ0::Union{Nothing,Matrix{ComplexF64}}=nothing)
    up = FiniteBath(HB_M, VB_M, β)
    dn = FiniteBath(HB_M, VB_M, β)
    if ρ0 !== nothing
        imp = TwoSpinImpurityModel(eps_d, eps_d, U, ρ0)
    elseif U == 0
        ρ1 = thermal2(eps_d, β)
        imp = TwoSpinImpurityModel(eps_d, eps_d, 0.0, Matrix{ComplexF64}(kron(ρ1, ρ1)))
    else
        ρt = thermal4(eps_d, U, β)
        imp = TwoSpinImpurityModel(eps_d, eps_d, U, ρt)
    end
    path = build_impurity_path(up, dn, N, T / N; order=order,
                               impurity_model=imp, cutoff=1e-12)
    return gtempo_greens_function_direct(path.factors..., path.lattice; flavor=1)
end

function compare(tag::String; order::Int, N::Int, bands::Int, U::Real,
                 ρ0::Union{Nothing,Matrix{ComplexF64}}=nothing)
    gt_n, lt_n = native_gtlt(order=order, N=N, bands=bands, U=U, ρ0=ρ0)
    gt_i, lt_i = ifgtempo_gtlt(order=order, N=N, U=U, ρ0=ρ0)
    ng = min(length(gt_n), length(gt_i))
    eg = maximum(abs.(gt_n[1:ng] .- gt_i[1:ng]))
    el = maximum(abs.(lt_n[1:ng] .- lt_i[1:ng]))
    ok = max(eg, el) < 1e-9
    @printf("%-34s : max|Δgt| = %.2e  max|Δlt| = %.2e  %s\n", tag, eg, el,
            ok ? "PASS" : "FAIL")
    return ok, gt_n, lt_n, gt_i, lt_i
end

function main()
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    all_ok = true
    println("═"^76)
    println("① 单带 Toulouse：native 紧凑 gt/lt vs IFGTEMPO 解耦双自旋 flavor-1")
    println("═"^76)
    for (order, N) in ((2, 3), (2, 4), (1, 3), (4, 1))
        ok, _, _, _, _ = compare("order $order N=$N (nl=$(order == 4 ? 3N : N))";
                                 order=order, N=N, bands=1, U=0.0)
        all_ok &= ok
    end
    println("═"^76)
    println("⑤ 双带 Anderson U=1.5：native bands=2 vs IFGTEMPO 直接同模型（热态 4×4）")
    println("═"^76)
    ρ0 = thermal4(EPS_D, UINT, BETA)
    for (order, N) in ((2, 3), (2, 4), (4, 1))
        ok, _, _, _, _ = compare("order $order N=$N U=1.5 (bands=2)";
                                 order=order, N=N, bands=2, U=UINT, ρ0=ρ0)
        all_ok &= ok
    end
    println("═"^76)
    println("附：双带 gt/lt 对 ED 连续演化的 Trotter 收敛（有效时刻）")
    println("═"^76)
    # 有效时刻：greater[k+1] ↔ (k+0.5)δt（k=0..N−1）；lesser[k+1] ↔ (k−0.5)δt（k=1..N）
    H, a, adag, H0 = singlemode_ed(μ=EPS_D, U=UINT, bathspecs=[(HB1, VB1^2), (HB1, VB1^2)],
                                   nbands=2)
    ρ = exp(-BETA .* H0); ρ = ρ / tr(ρ)
    UT = exp(-1im * T .* H)
    cH(op, t) = exp(-1im * t .* H)' * op * exp(-1im * t .* H)
    for (order, N) in ((2, 2), (2, 4))
        δt = T / N
        gt_n, lt_n = native_gtlt(order=order, N=N, bands=2, U=UINT)
        err_g, err_l = 0.0, 0.0
        for k in 0:(N - 1)
            τ = (k + 0.5) * δt
            ref = tr(cH(a, T) * cH(adag, τ) * ρ)     # ⟨a(T)a†(τ)⟩ 双 Heisenberg
            err_g = max(err_g, abs(gt_n[k + 1] - ref))
        end
        for k in 1:N
            τ = (k - 0.5) * δt
            ref = tr(cH(adag, τ) * cH(a, T) * ρ)     # ⟨a†(τ)a(T)⟩
            err_l = max(err_l, abs(lt_n[k + 1] - ref))
        end
        @printf("  order %d N=%d (δt=%.3f) : max|gt − gt_ed_eff| = %.3e  max|lt − lt_ed_eff| = %.3e\n",
                order, N, δt, err_g, err_l)
    end
    println("═"^76)
    println("总结：", all_ok ? "gt/lt 机器级对齐全绿（①单带 + ⑤双带）" : "存在 FAIL（见上）")
end

main()
