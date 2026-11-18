# ============================================================================
# bcsweyl/contact_extract_refined.jl — A₀ 的精化提取（二次拟合 + 大 N + 紧截断）
#
# 初版（contact_extract.jl）线性拟合被 O(δτ²) 项偏置 ~2%。本版：
#   g(δτ) = A₀·δτ + B·δτ² 二次拟合；等价地 g/δτ = A₀ + B·δτ 线性拟合取截距。
#   N ∈ {16,24,32,48,64}，IF 截断 D=256/rtol=1e-14（排除 IF 截断污染 1e-3 级）。
# 判据：cn*(Δ) = −4A₀/β 是否 ≡ 2α（Δ 无关，H1/Bogoliubov 正则性）。
# ============================================================================
include(joinpath(@__DIR__, "..", "test", "util.jl"))

using GTEMPO
using ImpurityModelBase
using LinearAlgebra
using Printf

const OUT = joinpath(@__DIR__, "results")
isdir(OUT) || mkpath(OUT)

const U, ϵ_ed, ω, α, β = 0.0, 0.0, 1.0, 0.5, 0.8
const Ns = (16, 24, 32, 48, 64)

function z_scheme(Nτ::Int, Δ; cn::Real)
    δτ = β / Nτ
    lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
    model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
    Is = bcs_hybriddynamics(lat, bath; trunc=Truncation(D=256, rtol=1.0e-14))
    K = sysdynamics_weyl2order(lat, model; coupling_norm=cn)
    real(integrate(K, Is))
end

function z_exact(Δ)
    H, Ab, Bb, H0 = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ)
    p1 = spin_half_matrices()
    n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
    I2 = one(n̂)
    Hbath = ω * (kron(n̂, I2) + kron(I2, n̂)) -
            (Δ * kron(JW * σ₊, σ₊) - conj(Δ) * kron(JW * σ₋, σ₋))
    real(tr(exp(-β * H))) / real(tr(exp(-β * Hbath)))
end

linfit(xs, ys) = begin
    n = length(xs); sx = sum(xs); sy = sum(ys)
    b = (n * sum(xs .* ys) - sx * sy) / (n * sum(abs2, xs) - sx^2)
    a = sy / n - b * sx / n
    a, b
end

io = open(joinpath(OUT, "contact_extract_refined.csv"), "w")
println(io, "delta,v2,u2,cn_star,minus2alpha_ratio,A0,B,max_resid")
println("="^96)
@printf("平凡系统 U=ϵ_d=0；期望检验：cn*(Δ) ≡ 2α = %.4f ？\n\n", 2α)
for Δ in (0.0, 0.3, 0.6, 1.0, 2.0)
    E = sqrt(ω^2 + abs2(Δ))
    u2 = (1 + ω / E) / 2
    v2 = (1 - ω / E) / 2
    Zed = z_exact(Δ)
    δτs = Float64[]; gs = Float64[]
    for N in Ns
        g = log(z_scheme(N, Δ; cn=0) / Zed)
        push!(δτs, β / N); push!(gs, g)
        @printf("  Δ=%.2f N=%2d  log gap=%+.8e  (g/δτ = %+.6f)\n", Δ, N, g, g / (β / N))
    end
    A0, B = linfit(δτs, gs ./ δτs)          # g/δτ = A₀ + B·δτ
    resid = maximum(abs, gs ./ δτs .- (A0 .+ B .* δτs))
    cn_star = -4 * A0 / β
    @printf(">>> Δ=%-4.2f cn* = %.6f  cn*/2α = %.5f  A₀ = %+.7f  B = %+.4f  resid = %.2e\n\n",
            Δ, cn_star, cn_star / (2α), A0, B, resid)
    println(io, "$Δ,$v2,$u2,$cn_star,$(cn_star / (2α)),$A0,$B,$resid")
end
close(io)
println("数据 → ", joinpath(OUT, "contact_extract_refined.csv"))
