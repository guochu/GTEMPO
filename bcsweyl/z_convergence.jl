# ============================================================================
# bcsweyl/z_convergence.jl — 裁决性基准：generic H_S 下 Z 的收敛阶
#
# contact_extract.jl（S=I 极限）给出所需 cn*(Δ)；本脚本在非平凡系统
# （U=1, ϵ_d≠0）上对候选 cn 验证 Z 的 log-log 收敛斜率：
#   order 1（参照）≈ 1；order 2 + 正确 cn ⟹ ≈ 2。
# 判据的 Δ 无关性：cn 候选若为真定义，须对所有 Δ 给出斜率 2。
#
# 运行：julia --project=<GTEMPO根> bcsweyl/z_convergence.jl [cn_mode]
#   cn_mode ∈ {sum_alpha, alpha, zero, plus_anom}（缺省 sum_alpha）
# ============================================================================
include(joinpath(@__DIR__, "..", "test", "util.jl"))   # bcs_ed / spin_half_matrices

using GTEMPO
using ImpurityModelBase
using LinearAlgebra
using Printf

const OUT = joinpath(@__DIR__, "results")
isdir(OUT) || mkpath(OUT)

# ── 非平凡系统 ──
const U, ϵ_ed, ω, α, β = 1.0, 0.5, 1.0, 0.5, 0.8
const DMAX = 192          # BCS IF 压缩上限（benchmarks 用默认 128；此处稍宽）
const CN_MODE = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :alpha

"候选 cn（K 侧 = 全带总耦合；两带 BCS 对 = 每带 α ⟹ 总 2α）：
 :sum_alpha = 2α（H1 预言：Bogoliubov 正则性，Δ 无关）；:alpha = α（半值对照）；
 :zero = 0（naive2 参照）；:plus_anom = 2α(u²+2v²)（反常道 2× 宽度计入候选）"
function candidate_cn(Δ)
    E = sqrt(ω^2 + abs2(Δ))
    u2 = (1 + ω / E) / 2
    v2 = (1 - ω / E) / 2
    CN_MODE === :sum_alpha && return 2α
    CN_MODE === :alpha && return α
    CN_MODE === :zero && return 0.0
    CN_MODE === :plus_anom && return 2α * (u2 + 2v2)   # 反常道按 2× 宽度计入
    error("unknown cn mode $CN_MODE")
end

function z_scheme(Nτ::Int, Δ; order::Int, cn::Real)
    δτ = β / Nτ
    lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
    model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
    bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
    Is = bcs_hybriddynamics(lat, bath; trunc=Truncation(D=DMAX, rtol=1e-12))
    K = order == 2 ?
        sysdynamics_weyl2order(lat, model; coupling_norm=cn) :
        sysdynamics(lat, model)
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

fit_slope(hs, errs) = begin
    x = log.(hs); y = log.(errs); n = length(hs)
    (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

const Ns = (8, 16, 24, 32, 48)
io = open(joinpath(OUT, "z_convergence_$(CN_MODE).csv"), "w")
println(io, "delta,order,N,dtau,Z,exact,abs_err")
println("="^92)
@printf("cn_mode = :%s\n", CN_MODE)
for Δ in (0.3, 0.6, 1.0)
    cn = candidate_cn(Δ)
    Zed = z_exact(Δ)
    @printf("Δ=%.1f  cn = %.6f (cn/α = %.4f)   Z_exact = %.8f\n", Δ, cn, cn / α, Zed)
    for order in (1, 2)
        errs = Tuple{Float64,Float64}[]
        for N in Ns
            Z = z_scheme(N, Δ; order=order, cn=cn)
            e = abs(Z - Zed)
            push!(errs, (β / N, e))
            println(io, "$Δ,$order,$N,$(β/N),$Z,$Zed,$e")
            @printf("  order=%d N=%2d  Z=%.8f  |Z−Zed| = %.3e\n", order, N, Z, e)
        end
        s = fit_slope(first.(errs), last.(errs))
        @printf("  >>> order=%d 斜率 = %.3f\n\n", order, s)
        println(io, "$Δ,$order,slope,,,,$s")
    end
end
close(io)
println("数据 → ", joinpath(OUT, "z_convergence_$(CN_MODE).csv"))
