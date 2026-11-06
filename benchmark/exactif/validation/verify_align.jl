# =============================================================================
# src/verify_align.jl —— 与 IFGTEMPO 的直接对齐验证 + quapi 路线缓存修复回归
#
# IFGTEMPO 对齐（U=0 解耦双自旋技巧）：两自旋杂质 U=0、两 bath 相同、初态
# ρth⊗ρth ⇒ 双自旋演化 = 单自旋演化 ⊗ 自身（精确恒等式）。IFGTEMPO GMPS
# 的 n_up（ρ[2,2]+ρ[4,4]，基序 |0⟩,|↑⟩,|↓⟩,|↑↓⟩）必须等于紧凑原生路线的
# n(T)（后者已与单自旋 Fock 链机器精度对齐）。
#
# quapi 回归：exactif_native_quapi 共享 _scheduled_gates!（缓存 bug 同受害者
# ——恒等门 gap 被上一 gap 的门张量替换），修复后地板应收紧。
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

const T_TOT = 0.4
const BETA = 1.3
const HB = ComplexF64[0.4;;]
const VV = ComplexF64[0.23 + 0.11im]
const EPS_D = 0.7

thermal_rho0(eps_d, β) = (H = ComplexF64[0 0; 0 eps_d]; exp(-β .* H) / tr(exp(-β .* H)))

"""紧凑原生路线的 n(T)（截断档生产用法）。"""
function compact_nt(order::Int, N::Int)
    model = ToulouseIM(ϵ_d=EPS_D)
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    pipe = NativeExact.exactif_native_compact(HB, VV, model, N, T_TOT / N; order=order,
                                              trunc=trunc, β=BETA, counter_c=-1)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k = lat.k
    mpsN = copy(mpsK)
    pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=1)]
    pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=1)]
    GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
    return real(GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z)
end

"""IFGTEMPO 两自旋 GMPS 的 n_up（U=0 解耦、相同 bath、ρth⊗ρth 初态）。"""
function ifgtempo_nt(order::Int, N::Int)
    up = FiniteBath(HB, VV, BETA)
    dn = FiniteBath(HB, VV, BETA)
    ρ1 = thermal_rho0(EPS_D, BETA)
    ρ0 = Matrix{ComplexF64}(kron(ρ1, ρ1))
    imp = TwoSpinImpurityModel(EPS_D, EPS_D, 0.0, ρ0)
    path = build_impurity_path(up, dn, N, T_TOT / N; order=order,
                               impurity_model=imp, cutoff=1e-12)
    ρ = gtempo_density_matrix_direct(path.factors...)
    n_up = real(ρ[2, 2] + ρ[4, 4])
    Z = real(tr(ρ))
    return n_up / Z
end

"""quapi 路线的 n(T)（修复后的回归；离散单模谱 = 有限模 bath 的等效谱）。
order 1 N=1 只有一个门/支（缓存 bug 惰性）——其与 Fock 的差 = 纯插入约定地板。
order 4 属已知 Δ 数据源限制（点基不承载跨 cell 门），不在回归范围。"""
function quapi_nt(order::Int, N::Int)
    bath = discretefermionicbath([0.4], [0.065]; β=BETA, μ=0.0)
    model = ToulouseIM(ϵ_d=EPS_D)
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    pipe = NativeExact.exactif_native_quapi(bath, model, N, T_TOT / N; order=order,
                                            trunc=trunc, β=BETA, counter_c=0)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k = lat.k
    mpsN = copy(mpsK)
    pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=1)]
    pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=1)]
    GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
    return real(GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z)
end

function main()
    if get(ARGS, 1, "") == "quapi"
        println("═"^76)
        println("quapi 路线回归（缓存修复后；离散单模谱）")
        println("═"^76)
        for (order, N) in ((1, 1), (2, 1), (2, 2))
            n_q = quapi_nt(order, N)
            up = FiniteBath(HB, VV, BETA)
            sched = trotter_schedule(order, N, T_TOT / N)
            ρ = fock_chain_single(up, EPS_D, sched, thermal_rho0(EPS_D, BETA))
            ref = real(ρ[2, 2]) / real(tr(ρ))
            @printf("order %d N=%d : quapi %.12f | Fock %.12f | Δ = %.3e\n",
                    order, N, n_q, ref, abs(n_q - ref))
        end
        return
    end
    println("═"^76)
    println("IFGTEMPO 直接对齐（两自旋 U=0 解耦恒等式；n_up[GMPS] vs n[紧凑原生]）")
    println("═"^76)
    for (order, N) in ((2, 1), (2, 2), (4, 1), (4, 2))
        n_ifg = ifgtempo_nt(order, N)
        n_nat = compact_nt(order, N)
        @printf("order %d N=%d : IFGTEMPO %.12f | 原生紧凑 %.12f | Δ = %.3e\n",
                order, N, n_ifg, n_nat, abs(n_ifg - n_nat))
    end
end

main()
