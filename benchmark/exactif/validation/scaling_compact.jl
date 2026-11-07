# =============================================================================
# src/scaling_compact.jl —— 紧凑消元路线的 Trotter 阶数标度（vs 精确连续演化）
#
# 参照 = 单模 bath 全 Hamiltonian（4×4）精确对角化的连续演化 ρ(T)；
# 误差 = |n_N − n_exact|，对 δt = T/N 做 log-log 最小二乘：
#   order 1 → 全局 p≈2（精确顶点幅重求和 [A,B]δt² 首项，同 IFGTEMPO 结论）
#   order 2 → p≈2；order 4 → p≈4（Suzuki 负 cell）
# 与 IFGTEMPO/twosector、FormB 的发表值（p = 2.00/2.00/4.05-4.16）对齐。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/src/scaling_compact.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(@__DIR__, "..", "src", "NativeExact.jl"))
using .NativeExact
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))

const T_TOT = 0.4
const BETA = 1.3
const HB = ComplexF64[0.4;;]
const VV = ComplexF64[0.23 + 0.11im]
const EPS_D = 0.7

thermal_rho0(eps_d, β) = (H = ComplexF64[0 0; 0 eps_d]; exp(-β .* H) / tr(exp(-β .* H)))

"""精确连续演化参照：单模 bath 全 H（4×4）精确对角化。"""
function exact_continuous()
    up = (bath_hamiltonian=HB, hybridization=VV, beta=BETA)
    m = 1
    nt = 1 + m
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt - 1)]
    hQ = NativeExact._coupled_h(VV, HB)             # 杂质能级 0（耦合单粒子 h）
    Hq = zeros(ComplexF64, dim, dim)
    for i in 1:(m + 1), j in 1:(m + 1)
        iszero(hQ[i, j]) || (Hq .+= hQ[i, j] .* (ops[i]' * ops[j]))
    end
    Him = EPS_D .* (ops[1]' * ops[1])
    Hb = zeros(ComplexF64, dim, dim)
    for i in 2:(m + 1), j in 2:(m + 1)
        iszero(HB[i - 1, j - 1]) || (Hb .+= HB[i - 1, j - 1] .* (ops[i]' * ops[j]))
    end
    Htot = Hq + Him   # hQ 已含 bath 能谱 + V 耦合（杂质能级 0），勿再加 Hb
    Zb = tr(exp(-BETA .* Hb)) / 2
    rb = exp(-BETA .* Hb) ./ Zb
    r0 = zeros(ComplexF64, dim, dim)
    ρ0 = thermal_rho0(EPS_D, BETA)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        r0[(b << 1 | a) + 1, (b << 1 | a2) + 1] = ρ0[a + 1, a2 + 1]
    end
    U = exp(-1im * T_TOT .* Htot)
    rtot = U * (r0 * rb) * adjoint(U)
    ρ = zeros(ComplexF64, 2, 2)
    for b in 0:((1 << m) - 1), a in 0:1, a2 in 0:1
        ρ[a + 1, a2 + 1] += rtot[(b << 1 | a) + 1, (b << 1 | a2) + 1]
    end
    return real(ρ[2, 2]) / real(tr(ρ))
end

_slope(hs, es) = (sum((hs .- sum(hs) / length(hs)) .* log.(es)) /
                  sum(abs2, hs .- sum(hs) / length(hs)))   # ln err vs ln δt → +p

function main()
    n_ex = exact_continuous()
    println("═"^76)
    println("紧凑路线 Trotter 标度 vs 精确连续演化（T=$(T_TOT), β=$(BETA), εd=$(EPS_D)）")
    @printf("连续精确 n(T) = %.12f\n", n_ex)
    println("═"^76)
    trunc = GTEMPO.truncdimcutoff(D=512, ϵ=1e-12)
    Ns = [1, 2, 4, 8]
    results = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for order in (1, 2, 4)
        pts = Tuple{Float64,Float64}[]
        for N in Ns
            model = ToulouseIM(ϵ_d=EPS_D)
            pipe = NativeExact.exactif_native_compact(HB, VV, model, N, T_TOT / N;
                                                      order=order, trunc=trunc, β=BETA,
                                                      counter_c=-1)
            (; lat, mpsK, mpsIs, norm) = pipe
            Z = GTEMPO.integrate(lat, mpsK, mpsIs...) * norm
            k = lat.k
            mpsN = copy(mpsK)
            pos_a = lat[GrassmannContourIndex(k, conj=false, branch=:+, band=1)]
            pos_ab = lat[GrassmannContourIndex(k, conj=true, branch=:-, band=1)]
            GTEMPO.apply!(GTEMPO.GTerm(pos_a, pos_ab, coeff=1.0), mpsN)
            nt = real(GTEMPO.integrate(lat, mpsN, mpsIs...) * norm / Z)
            err = abs(nt - n_ex)
            δt = T_TOT / N
            push!(pts, (δt, err))
            @printf("order %d N=%2d (δt=%.4f) : n = %.12f | err = %.3e\n",
                    order, N, δt, nt, err)
        end
        results[order] = pts
    end
    println("═"^76)
    println("最小二乘阶数（err ∝ δtᵖ）：")
    for order in (1, 2, 4)
        pts = results[order]
        p = _slope([log(δt) for (δt, e) in pts], [e for (δt, e) in pts])
        @printf("order %d : p = %.2f\n", order, p)
    end
end

main()
