# =============================================================================
# src/verify_align_u15.jl —— IFGTEMPO 直接同模型对齐（U=1.5 Anderson，bands=2）
#
# 与 verify_align.jl 的 U=0 解耦技巧不同，这里 IFGTEMPO 两自旋路线直接跑
# 同一个交互模型：TwoSpinImpurityModel(0.7, 0.7, U=1.5) + 相同单模 bath ×2 +
# 相同热态 ρ₀（4×4，含 U 关联）。可观测量两个（独立信息）：
#   n_up  = (ρ₂₂+ρ₄₄)/tr   —— band-1 折返腿 GTerm 投影
#   ⟨n₁n₂⟩ = ρ₄₄/tr        —— band-1+band-2 两个折返腿 GTerm 的乘积投影
# 先用 Fock 离散链（16 维直算）机器级校验原生读出约定，再对 IFGTEMPO。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/verify_align_u15.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(@__DIR__, "..", "src", "NativeExact.jl"))
using .NativeExact
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))

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

# ---- Fock 离散链（16 维；mode 序 imp₁=0, imp₂=1, bath₁=2, bath₂=3）----
function annih_op(dim, mode)
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    return op
end

function chain_n(n1n2::Bool, order::Int, N::Int)
    ops = [annih_op(16, k) for k in 0:3]
    n = [ops[k]' * ops[k] for k in 1:4]
    Hq = zeros(ComplexF64, 16, 16)
    for (imp, bath) in ((1, 3), (2, 4))
        Hq .+= HB1 .* n[bath]
        Hq .+= VB1 .* (ops[imp]' * ops[bath])
        Hq .+= VB1 .* (ops[bath]' * ops[imp])
    end
    Him = EPS_D .* (n[1] + n[2]) .+ UINT .* (n[1] * n[2])
    rb1 = exp(-BETA .* HB1 .* ComplexF64[0 0; 0 1]); rb1 = rb1 / tr(rb1)
    rb = zeros(ComplexF64, 16, 16)
    for b1 in 0:1, b2 in 0:1, b1p in 0:1, b2p in 0:1
        b = (b2 << 1 | b1); bp = (b2p << 1 | b1p)
        for a in 0:3
            rb[(b << 2 | a) + 1, (bp << 2 | a) + 1] = rb1[b1 + 1, b1p + 1] * rb1[b2 + 1, b2p + 1]
        end
    end
    ρ0 = thermal4(EPS_D, UINT, BETA)
    r0 = zeros(ComplexF64, 16, 16)
    for b in 0:3, a in 0:3, a2 in 0:3
        r0[(b << 2 | a) + 1, (b << 2 | a2) + 1] = ρ0[a + 1, a2 + 1]
    end
    sched = trotter_schedule(order, N, T / N)
    U = chain_operator(Hq, Him, collect(Float64.(sched.cells)),
                       collect(Float64.(sched.links)))
    rt = U * (r0 * rb) * adjoint(U)
    obs = n1n2 ? (n[1] * n[2]) : n[1]
    return real(tr(rt * obs)) / real(tr(rt))
end

# ---- 原生紧凑路线（bands=2）：n_up = band-1 折返腿 GTerm；⟨n₁n₂⟩ = 双 GTerm ----
function native_obs(n1n2::Bool, order::Int, N::Int)
    model = GTEMPO.AndersonIM(U=UINT, ϵ_d=EPS_D)
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    pipe = NativeExact.exactif_native_compact(HB_M, VB_M, model, N, T / N;
                                              order=order, trunc=trunc, β=BETA,
                                              bands=2, counter_c=-1)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k = lat.k
    mpsN = copy(mpsK)
    for band in (n1n2 ? (1, 2) : (1,))
        pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=band)]
        pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=band)]
        GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
    end
    return real(GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z)
end

# ---- IFGTEMPO 直接同模型（U=1.5，热态 4×4 初态，相同 bath ×2）----
function ifgtempo_obs(n1n2::Bool, order::Int, N::Int)
    up = FiniteBath(HB_M, VB_M, BETA)
    dn = FiniteBath(HB_M, VB_M, BETA)
    imp = TwoSpinImpurityModel(EPS_D, EPS_D, UINT, thermal4(EPS_D, UINT, BETA))
    path = build_impurity_path(up, dn, N, T / N; order=order,
                               impurity_model=imp, cutoff=1e-12)
    ρ = gtempo_density_matrix_direct(path.factors...)
    trρ = real(tr(ρ))
    return n1n2 ? real(ρ[4, 4]) / trρ : real(ρ[2, 2] + ρ[4, 4]) / trρ
end

function main()
    println("═"^78)
    println("① 原生读出约定机器级校验（紧凑管线 vs Fock 离散链；U=1.5 热态）")
    println("═"^78)
    for (order, N) in ((2, 1), (2, 2), (4, 1), (4, 2)), (tag, n1n2) in (("n_up", false), ("n₁n₂", true))
        ref = chain_n(n1n2, order, N)
        nat = native_obs(n1n2, order, N)
        @printf("  order %d N=%d %-5s : pipeline %.13f | chain %.13f | Δ = %.2e\n",
                order, N, tag, nat, ref, abs(nat - ref))
    end
    println("═"^78)
    println("② IFGTEMPO 直接同模型对齐（U=1.5，双可观测量）")
    println("═"^78)
    all_ok = true
    for (order, N) in ((2, 1), (2, 2), (4, 1), (4, 2)), (tag, n1n2) in (("n_up", false), ("n₁n₂", true))
        nat = native_obs(n1n2, order, N)
        ifg = ifgtempo_obs(n1n2, order, N)
        ok = abs(nat - ifg) < 1e-10
        all_ok &= ok
        @printf("  order %d N=%d %-5s : IFGTEMPO %.13f | 原生紧凑 %.13f | Δ = %.2e  %s\n",
                order, N, tag, ifg, nat, abs(nat - ifg), ok ? "PASS" : "FAIL")
    end
    println("═"^78)
    println("总结：", all_ok ? "U=1.5 直接同模型对齐全绿（双可观测量）" : "存在 FAIL")
end

main()
