# =============================================================================
# src/verify_compact.jl —— 紧凑消元路线（exactif_native_compact）的对齐验证
#
# 判据（ntT 同款折返腿 GTerm 投影读出，Z = integrate·ν·Pf）：
#   ① 离散链对齐：n(T) vs Fock 离散 Trotter 链（同 schedule、同热态 ρ₀）；
#   ② Z 证书：εd = 0 时 Z = integrate·ν·Pf ≈ 1；
#   ③ V = 0 证书：n(T) = n₀（热态期望）对任意 order/N。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/verify_compact.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: TruncationScheme, NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(@__DIR__, "NativeExact.jl"))
using .NativeExact

include(joinpath(GT_ROOT, "benchmark", "exactif", "form_b", "formb_dense.jl"))

const T_TOT = 0.4
const BETA = 1.3
const HB = ComplexF64[0.4;;]
const VV = ComplexF64[0.23 + 0.11im]
const EPS_D = 0.7

"""热态 ρ₀ = e^{-βH}/tr（与 systhermalstate! 同一语义）。"""
function thermal_rho0(eps_d::Real, β::Real)
    H = ComplexF64[0 0; 0 eps_d]
    R = exp(-β .* H)
    return R / tr(R)
end

"""紧凑管线的 n(T)（ntT 同款：折返点自由腿对 (a⁺_k, ā⁻_k) 的 2 腿 GTerm 投影）。

Z = integrate·norm（norm = ν·Pf，管线内置）。"""
function compact_nt(; order::Int, N::Int, eps_d::Real=EPS_D, V::Vector{ComplexF64}=VV,
                    counter_c::Int=-1, trunc::TruncationScheme=NoTruncation(),
                    β::Real=BETA)
    model = ToulouseIM(ϵ_d=eps_d)
    pipe = NativeExact.exactif_native_compact(HB, V, model, N, T_TOT / N; order,
                                              trunc, β, counter_c)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k = lat.k
    mpsN = copy(mpsK)
    pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=1)]
    pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=1)]
    GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
    nt = real(GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z)
    return (Z = Z, nt = nt, k = k, nl = length(pipe.cells))
end

"""Fock 离散链参照（同 schedule、同热态 ρ₀）。"""
function fock_ref(order::Int, N::Int, eps_d::Real)
    up = FiniteBath(HB, VV, BETA)
    sched = trotter_schedule(order, N, T_TOT / N)
    ρ = fock_chain_single(up, eps_d, sched, thermal_rho0(eps_d, BETA))
    return real(ρ[2, 2]) / real(tr(ρ))
end

function main(; tier::Symbol=:exact)
    # 精确档（NoTruncation）：orders 1/2 × N≤3 + order 4 × N=1（消除核键 4^nl 增长，
    # order 4 N≥2 必须走截断档）；截断档（D=512, ϵ=1e-12，twosector 生产用法）：
    # order 4 × N=2,3。
    cases = tier === :exact ?
            Tuple{Int,Int}[(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (4, 1)] :
            Tuple{Int,Int}[(1, 4), (2, 4), (4, 2), (4, 3)]
    trunc = tier === :exact ? NoTruncation() :
            GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    println("═"^76)
    println("紧凑消元路线 n(T) vs Fock 离散链（$tier 档；T=$(T_TOT), β=$(BETA), εd=$(EPS_D)）")
    println("═"^76)
    all_ok = true
    for (order, N) in cases
        ref = fock_ref(order, N, EPS_D)
        try
            out = compact_nt(order=order, N=N, trunc=trunc)
            ok = abs(out.nt - ref) < 1e-9
            all_ok &= ok
            @printf("order %d N=%d (nl=%d, k=%d) : n = %.10f | Fock %.10f | Δn = %.3e  %s\n",
                    order, N, out.nl, out.k, out.nt, ref, abs(out.nt - ref),
                    ok ? "PASS" : "FAIL")
        catch e
            all_ok = false
            @printf("order %d N=%d : 管线异常：%s\n", order, N,
                    sprint(showerror, e)[1:80])
        end
    end

    println("═"^76)
    println("Z 证书（εd = 0：Z = integrate·ν·Pf 应为 1）")
    for (order, N) in (tier === :exact ?
                       ((1, 1), (2, 1), (2, 2), (4, 1)) : ((4, 2),))
        try
            out = compact_nt(order=order, N=N, eps_d=0.0, trunc=trunc)
            @printf("order %d N=%d : Z−1 = %.3e\n", order, N, abs(out.Z - 1))
        catch e
            @printf("order %d N=%d : 管线异常：%s\n", order, N,
                    sprint(showerror, e)[1:80])
        end
    end

    println("═"^76)
    println("V = 0 证书（n(T) = n₀ 热态期望，任意 order/N）")
    n0 = real(thermal_rho0(EPS_D, BETA)[2, 2])
    @printf("n₀ = %.12f\n", n0)
    for (order, N) in (tier === :exact ? ((1, 1), (2, 2), (4, 1)) : ((4, 2),))
        try
            out = compact_nt(order=order, N=N, V=ComplexF64[0.0 + 0.0im], trunc=trunc)
            @printf("order %d N=%d : n = %.14f | n₀ %.14f | Δ = %.3e\n",
                    order, N, out.nt, n0, abs(out.nt - n0))
        catch e
            @printf("order %d N=%d : 管线异常：%s\n", order, N,
                    sprint(showerror, e)[1:80])
        end
    end
    println("═"^76)
    println("总结（$(tier) 档）：", all_ok ? "离散链对齐全绿" : "存在 FAIL（见上）")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(tier = contains(string(get(ARGS, 1, "")), "trunc") ? :trunc : :exact)
end
