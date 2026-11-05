# =============================================================================
# src/verify_native.jl —— NativeExact（GTEMPO 原生消元式管线）的对齐验证
#
# 判据：
#   ① Z 证书：ε_d = 0 时 Z = integrate(lat, mpsK, mpsIs…)·norm == tr(ρ₀) = 1；
#   ② Fock 判决：ε_d ≠ 0 时 ⟨n(T)⟩ 与 greater/lesser 端点 vs Fock 离散链；
#   ③（后续）两带 vs twosector / FormB 参照、截断档、Trotter 标度。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/verify_native.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: TruncationScheme, truncdimcutoff, NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(@__DIR__, "NativeExact.jl"))
using .NativeExact

# Fock 参照（formb_dense 的离散链机器；参照端允许 IFGTEMPO）
include(joinpath(GT_ROOT, "benchmark", "exactif", "form_b", "formb_dense.jl"))

const T_TOT = 0.4
const BETA = 1.3
const HB = ComplexF64[0.4;;]
const VV = ComplexF64[0.23 + 0.11im]

_gtempo_model(eps_d) = ToulouseIM(ϵ_d=eps_d)

function run_case(; order::Int, N::Int, eps_d::Real, counter_c::Int,
                  trace_style::Symbol, trunc::TruncationScheme)
    sched = native_schedule(order, N, T_TOT / N)
    model = _gtempo_model(eps_d)
    pipe = exactif_native_pipeline(((hb = HB, v = VV),), model, N, T_TOT / N;
                                   order, trunc, β = BETA, bands = 1,
                                   trace_style, counter_c)
    (; lat, mpsK, mpsIs, norm) = pipe
    Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
    k = lat.k
    nt = GTEMPO.occupation(lat, k, mpsK, mpsIs...; branch = :+, Z = Z)
    gt = GTEMPO.greater(lat, k, 1, mpsK, mpsIs...; Z = Z)
    lt = GTEMPO.lesser(lat, 1, k, mpsK, mpsIs...; Z = Z)
    return (Z = Z, nt = real(nt), gt = gt, lt = lt)
end

function main(; Nmax::Int=2)
    println("═"^78)
    println("NativeExact 首轮标定（单自旋，T=$(T_TOT)，β=$(BETA)，单模 bath）")
    println("═"^78)
    trunc = NoTruncation()

    # ---- Fock 参照 ----
    up = FiniteBath(HB, VV, BETA)
    imp = (eps = 0.7, rho0 = ComplexF64[1 0; 0 0])
    refs = Dict{Tuple{Int,Int},Any}()
    for order in (1, 2), N in 1:Nmax
        sched = trotter_schedule(order, N, T_TOT / N)
        refs[(order, N)] = fock_chain_single(up, 0.7, sched, imp.rho0)
    end

    best = Dict{Tuple{Int,Int},NamedTuple}()
    for counter_c in (1, -1), trace_style in (:edge, :flip)
        ok_all = true
        for order in (1, 2), N in 1:Nmax
            # Z 证书（εd = 0）
            out0 = run_case(order = order, N = N, eps_d = 0.0,
                            counter_c = counter_c, trace_style = trace_style,
                            trunc = trunc)
            errz = abs(out0.Z - 1)
            # Fock 判决（εd ≠ 0）
            out = run_case(order = order, N = N, eps_d = 0.7,
                           counter_c = counter_c, trace_style = trace_style,
                           trunc = trunc)
            ref = refs[(order, N)]
            n_ref = real(ref[2, 2]) / real(tr(ref))
            errn = abs(out.nt - n_ref)
            score = errz + errn
            if !haskey(best, (order, N)) || score < best[(order, N)].score
                best[(order, N)] = (counter_c = counter_c,
                                    trace_style = trace_style, score = score,
                                    errz = errz, errn = errn, nt = out.nt,
                                    n_ref = n_ref)
            end
            @printf("  order %d N=%d cc=%+d %-5s : Z−1 = %9.2e  ⟨n⟩ = %.8f (Fock %.8f, err %.2e)\n",
                    order, N, counter_c, trace_style, errz, out.nt, n_ref, errn)
            ok_all &= (errz < 1e-9 && errn < 1e-9)
        end
        ok_all && println("  ↑ cc=$counter_c $trace_style 全部通过！")
    end
    println("═"^78)
    for order in (1, 2), N in 1:Nmax
        b = best[(order, N)]
        @printf("最优 order %d N=%d : cc=%+d %-5s  Z−1 = %.2e  n-err = %.2e\n",
                order, N, b.counter_c, b.trace_style, b.errz, b.errn)
    end
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
