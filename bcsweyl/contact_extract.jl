# ============================================================================
# bcsweyl/contact_extract.jl — BCS 浴 Weyl 标量系数的接触缺口提取
#
# 原理（推导见 README.md §2）：
#   order=2 K 的 Weyl 标量 exp(cn/4·Σ δτ²) 是全局 O(δτ) 相对因子：
#     log Z(cn) = log Z_naive2(cn=0) + cn·β·δτ/4 + O(δτ²)   （τ 环 N 格，δτ=β/N）
#   取平凡系统极限 U = ϵ_d = 0 ⟹ 传播子 S = I，order1/2 之差纯为
#   {q=1/2 权重平移, W₄ I/2, Weyl 标量}——接触机器被完全隔离。
#   naive2 的 O(δτ) 对数缺口系数 A₀ 直接给出所需 cn*：
#     log(Z_naive2/Z_exact) = A₀·δτ + O(δτ²)  ⟹  cn* = −4·A₀/β
#
# 机器验证：Δ=0（BCS bath 退化为两条独立 normal 带）必须给出 cn* = α。
# 测量对象：BCS bath 的 cn*(Δ)——若与 Δ 无关且 = α，则
#   coupling_norm(::AbstractBCSBath) = coupling_norm(裸谱)（Bogoliubov 正则性）；
# 若随 Δ 变化，则反常道进入标量，定义须含 Δ 依赖项。
#
# 运行：julia --project=<GTEMPO根> bcsweyl/contact_extract.jl
# ============================================================================
include(joinpath(@__DIR__, "..", "test", "util.jl"))   # bcs_ed / spin_half_matrices

using GTEMPO
using ImpurityModelBase
using LinearAlgebra
using Printf

const OUT = joinpath(@__DIR__, "results")
isdir(OUT) || mkpath(OUT)

# ── 平凡系统参数：H_imp = 0 ⟹ S = I（order1/2 只差接触机器）──
const U, ϵ_ed, ω, α, β = 0.0, 0.0, 1.0, 0.5, 0.8

"Z(K, IF) at order 2 with coupling_norm cn (trivial system, imag ring)"
function z_scheme(Nτ::Int, Δ; cn::Real)
    δτ = β / Nτ
    lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
    model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
    Is = bcs_hybriddynamics(lat, bath)
    K = sysdynamics_weyl2order(lat, model; coupling_norm=cn)
    real(integrate(K, Is))
end

"ED 参照（bcs_ed 哈密顿量 / 浴自身配对配分——test/bcs 同款归一化）"
function z_exact(Δ)
    H, Ab, Bb, H0 = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ)
    p1 = spin_half_matrices()
    n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
    I2 = one(n̂)
    Hbath = ω * (kron(n̂, I2) + kron(I2, n̂)) -
            (Δ * kron(JW * σ₊, σ₊) - conj(Δ) * kron(JW * σ₋, σ₋))
    real(tr(exp(-β * H))) / real(tr(exp(-β * Hbath)))
end

"线性拟合 y = a·x + b 的 a"
slope_fit(xs, ys) = begin
    n = length(xs); sx = sum(xs); sy = sum(ys)
    (n * sum(xs .* ys) - sx * sy) / (n * sum(abs2, xs) - sx^2)
end

"""
    extract_cn_star(Δ; Ns) -> (cn*, A₀, resid)

naive2 (cn=0) 的 log 缺口 g(δτ) = A₀·δτ + O(δτ²) 线性拟合 ⟹ cn* = −4·A₀/β。
"""
function extract_cn_star(Δ; Ns=(8, 16, 24, 32, 40))
    Zed = z_exact(Δ)
    δτs = Float64[]; gs = Float64[]
    for N in Ns
        g = log(z_scheme(N, Δ; cn=0) / Zed)
        push!(δτs, β / N); push!(gs, g)
        @printf("    Δ=%.2f N=%2d  Z_naive2=%.8f  Z_exact=%.8f  log gap=%+.6e\n",
                Δ, N, z_scheme(N, Δ; cn=0), Zed, g)
    end
    A0 = slope_fit(δτs, gs)
    cn_star = -4 * A0 / β
    resid = maximum(abs, gs .- A0 .* δτs)
    return cn_star, A0, resid
end

# ── 主流程 ─────────────────────────────────────────────────────────────────
io = open(joinpath(OUT, "contact_extract.csv"), "w")
println(io, "delta,cn_star,A0,max_resid")
println("="^88)
println("机器验证（Δ=0 ⟹ 两条独立 normal 带，期望 cn* = α = $α）")
for Δ in (0.0, 0.3, 0.6, 1.0, 2.0)
    cn_star, A0, resid = extract_cn_star(Δ)
    tag = Δ == 0.0 ? "  (normal 验证)" : ""
    @printf(">>> Δ=%-4.2f cn* = %.6f  (cn*/α = %.4f)  A₀ = %+.6e  max_resid = %.2e%s\n",
            Δ, cn_star, cn_star / α, A0, resid, tag)
    println(io, "$Δ,$cn_star,$A0,$resid")
end
close(io)
println("="^88)
println("数据 → ", joinpath(OUT, "contact_extract.csv"))
