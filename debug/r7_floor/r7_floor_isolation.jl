# =============================================================================
# R7 排查实验：Weyl2 order=2 ~1e-5 δt 无关地板的隔离
# 四臂对照（同一 harness：DiracDelta 单模 ω=1 α=0.5 β=1，Anderson U=1 ϵ_d=0.8，
# T=0.4，greater/lesser vs ED，时刻对齐切片）：
#   A: bcsbath Δ=0.6, bands=2 :z2   （原基准——地板复现）
#   B: bcsbath Δ=0.0, bands=2 :z2   （纯 normal 2 带——反常通道排除）
#   C: normal bath,   bands=2       （normal 2 带同框架——BCS 路径排除）
#   D: normal bath,   bands=1       （九组基准类——带数对照）
# 判据：地板出现在哪些臂 ⟹ 根因面（带数 vs BCS 路径 vs 框架）。
# =============================================================================
using GTEMPO, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "test", "util.jl"))
include(joinpath(@__DIR__, "..", "test", "normalbath", "util.jl"))

const U = 1.0; const ϵ_ed = 0.8; const ω = 1.0; const α = 0.5
const β = 1.0; const T = 0.4
const NTS = (4, 8, 16, 32)
const NMAX = 32

p1 = spin_half_matrices()
n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
I2 = one(n̂)
n_imp = kron(n̂, I2) + kron(I2, n̂)
Himpbare = -ϵ_ed * n_imp + U * kron(n̂, n̂)
I4 = Matrix{ComplexF64}(I, 4, 4)
n_bath = kron(I2, n̂) + kron(n̂, I2)

"跑一个臂：返回 Dict((Nt) => (eg, el))"
function run_arm(label; bcs::Bool, bands::Int, Δ=0.6, sym=Z2Symmetry())
    ts = (1:NMAX) .* (T / NMAX)
    model = bands == 1 ? ToulouseIM(ϵ_d=-ϵ_ed) : AndersonIM(U=U, ϵ_d=-ϵ_ed)
    ρ0 = thermal_impurity_state(model, β)
    if bcs
        # R7 修正：ED 初态 = 杂质热态 ⊗ **配对浴**热态（bcs_ed 返回的 H0 含浴配对项
        # ——TN 的 kernel 编码配对浴；未配对 H0dec 曾制造 ~1e-5 δt⁰ 假地板）
        H, Ab, Bb, H0bcs = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ)
        gt_ed, lt_ed = greater_lesser_ed(H, Ab, Bb, H0bcs, ts, β)
        bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
        cn = coupling_norm(bath) * bands
    else
        H, a, adag, H0 = singlemode_ed(μ=-ϵ_ed, U=U; bathspecs=[(ω, α) for _ in 1:bands])
        gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
        bath = discretefermionicbath([ω], [α]; β=β)
        cn = coupling_norm(bath) * bands
    end
    errs = Dict{Int,Tuple{Float64,Float64}}()
    for Nt in NTS
        δt = T / Nt
        kw = bands == 1 ? () : (bands=bands, symmetry=sym)
        lat = GrassmannLattice(N=Nt, δt=δt, contour=:real; kw...)
        kernel = correlationfunction(bath, lat)
        Is, _ = hybriddynamics(lat, kernel, PartialIF())
        Is = Is isa GMPS ? (Is,) : Is
        K = sysdynamics(lat, model; order=2, coupling_norm=cn, initial_state=ρ0)
        Z = integrate(K, Is...)
        stride = NMAX ÷ Nt
        idx = (1:Nt) .* stride
        gt = [greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2,
                      coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
        lt = [lesser(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2,
                     coupling_norm=cn, initial_state=ρ0) for i in 1:Nt]
        errs[Nt] = (maximum(abs.(gt .- gt_ed[idx])), maximum(abs.(lt .- lt_ed[idx])))
    end
    println("\n== 臂 $label（bcs=$bcs, bands=$bands, Δ=$Δ, cn=$cn）==")
    println("  Nt    maxabs_gt      maxabs_lt")
    for Nt in NTS
        @printf("  %3d   %.3e    %.3e\n", Nt, errs[Nt][1], errs[Nt][2])
    end
    for (idx, lbl) in ((1, "gt"), (2, "lt"))
        x = [log(T / Nt) for Nt in NTS]
        y = [log(max(errs[Nt][idx], 1e-16)) for Nt in NTS]
        n = length(x)
        slope = sum((x[i] - sum(x)/n) * (y[i] - sum(y)/n) for i in 1:n) /
                sum((x[i] - sum(x)/n)^2 for i in 1:n)
        @printf("  slope(%s) = %.3f\n", lbl, slope)
    end
    return errs
end

run_arm("A: BCS Δ=0.6"; bcs=true, bands=2, Δ=0.6)
run_arm("B: BCS Δ=0"; bcs=true, bands=2, Δ=0.0)
run_arm("C: normal bands=2"; bcs=false, bands=2)
run_arm("D: normal bands=1 (U=0 Toulouse)"; bcs=false, bands=1)
