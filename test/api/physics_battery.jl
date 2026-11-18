# Phase 5c 物理电池：M=2 双带 Anderson × {U=0, U=0.8} × {real, imag} 轮廓
# 对照 8 维精确对角化（ED）。四重验证合一：
#   ① 带因子化（U=0 时 M=2 ≡ M=1，机器精度）；
#   ② ED 收敛（δt/δτ 减半误差按二阶收缩）；
#   ③ 跨对称性一致（:z2/:u1/:z2z2/:u1u1 四形逐位相等）；
#   ④ 浴落带正确性（band-2 无浴时 ⟨n₂⟩ = bathless 常数精确）。
#
# ED 约定（⚠️ 三次踩坑后定谳，勿改）：浴模在**绝对 ω，无 μ 移**——
#   H_ED = μΣn + U n₁n₂ + ω·m + √α(b†c₁ + c₁†b)，band 1 = 浴耦合带；
#   实时：ρ₀(杂质热态)⊗ρ_bath(浴热态) 演化；虚时：全 8 维热态。
# 对照设定：只有 band-1 带 IF、band-2 真空链（U≠0 时双带浴的 U 反馈使
# 8 维 ED 不适用——须 16 维双浴 ED 或真空对照）。
#
# 用法：julia --project=. test/api/physics_battery.jl

using Test, Printf
using GTEMPO
using ImpurityModelBase
using LinearAlgebra: I, eigen, tr, Diagonal, Hermitian, kron

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
const bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
const BATHLESS = exp(-β * μ) / (1 + exp(-β * μ))

# ---------------------------------------------------------------- ED 参照 --
function _annih(m::Int, modes::Int)
    d = zeros(ComplexF64, 2^modes, 2^modes)
    for st in 0:(2^modes)-1
        bit = 1 << (m - 1)
        (st & bit) != 0 || continue
        d[st ⊻ bit + 1, st + 1] = iseven(count_ones(st & (bit - 1))) ? 1.0 : -1.0
    end
    d
end

"8 维 ED（band 1 耦浴）。kind = :thermal（全热态 ⟨n₁⟩）/ :real（演化 ⟨n₁(t)⟩）"
function ed_reference(U; kind, dt=0.0, Nt=0)
    I2 = Matrix{ComplexF64}(I, 2, 2)
    c1 = kron(_annih(1, 2), I2)
    c2 = kron(_annih(2, 2), I2)
    bm = kron(Matrix{ComplexF64}(I, 4, 4), _annih(1, 1))
    n1, n2 = c1' * c1, c2' * c2
    nm = bm' * bm
    H = μ * (n1 + n2) + U * (n1 * n2) + ω * nm + sqrt(α) * (bm' * c1 + c1' * bm)
    ev, vec = eigen(Hermitian(Matrix(H)))
    w = exp.(-β .* (ev .- minimum(ev))); w ./= sum(w)
    ρ = vec * Diagonal(ComplexF64.(w)) * vec'
    kind === :thermal && return real(tr(ρ * n1))
    p = [1.0, exp(-β * μ), exp(-β * μ), exp(-β * (2μ + U))]
    ρ0 = Diagonal(ComplexF64.(p)); ρ0 /= tr(ρ0)
    evb, vecb = eigen(Hermitian(Matrix(ω * [0.0 0.0; 0.0 1.0])))
    wb = exp.(-β .* (evb .- minimum(evb))); wb ./= sum(wb)
    ρb = vecb * Diagonal(ComplexF64.(wb)) * vecb'
    ρ = kron(ρ0, ρb)
    Ue = exp(-1.0im * dt * H)
    for _ in 1:Nt; ρ = Ue * ρ * Ue'; end
    return real(tr(ρ * n1))
end

# ---------------------------------------------------------------- 电池 --
@testset "Phase 5c 物理电池：M=2 vs ED" begin
    for U in (0.0, 0.8)
        # ---------- 虚时：全热态 ----------
        ed = ed_reference(U; kind=:thermal)
        last_err = Inf
        for (N, dτ) in ((4, 0.25), (8, 0.125))
            lat = GrassmannLattice(N=N, δτ=dτ, bands=2, contour=:imag)
            cn = α * 2
            kernel = correlationfunction(bath, lat; band=1)
            IF1, _ = hybriddynamics(lat, kernel, PartialIF())
            vac = vacuumstate(lat)
            K = sysdynamics(lat, AndersonIM(U=U, ϵ_d=μ); order=2, coupling_norm=cn)
            zs = integrate_scaled(K, IF1, vac)
            Z = exp(zs.logabs) * zs.phase
            occ = occupation(lat, AndersonIM(U=U, ϵ_d=μ), nothing, IF1, vac; j=N ÷ 2, branch=:τ,
                             band=1, order=2, coupling_norm=cn, Z=Z)
            err = abs(real(occ) - ed)
            @printf("  imag U=%.1f N=%-2d: ⟨n₁⟩=%.7f  ED=%.7f  err=%.2e\n", U, N, real(occ), ed, err)
            @test err < 1e-3                     # 分辨率锚（δτ=0.25 二阶 ≈1e-3 档）
            N == 8 && (last_err = err)
            # ④ band-2 无浴 → bathless 常数（τ-不变）
            occ2 = occupation(lat, AndersonIM(U=U, ϵ_d=μ), nothing, IF1, vac; j=N ÷ 2, branch=:τ,
                              band=2, order=2, coupling_norm=cn, Z=Z)
            U == 0 && @test abs(real(occ2) - BATHLESS) < 1e-10
        end
        @test last_err < 4e-4                    # δτ=0.125 收敛锚
        # ---------- 跨对称性一致（:z2 基准 vs 其余三形逐位相等）----------
        N = 4
        lat0 = GrassmannLattice(N=N, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
        cn = α * 2
        kernel0 = correlationfunction(bath, lat0; band=1)
        IF0, _ = hybriddynamics(lat0, kernel0, PartialIF())
        vac0 = vacuumstate(lat0)
        K0 = sysdynamics(lat0, AndersonIM(U=U, ϵ_d=μ); order=2, coupling_norm=cn)
        zs0 = integrate_scaled(K0, IF0, vac0)
        ref_Z = exp(zs0.logabs) * zs0.phase
        ref_occ = occupation(lat0, AndersonIM(U=U, ϵ_d=μ), nothing, IF0, vac0; j=2, branch=:τ,
                             band=1, order=2, coupling_norm=cn, Z=ref_Z)
        for sym in (U1Symmetry(), Z2Z2Symmetry(), U1U1Symmetry())
            lat = GrassmannLattice(N=N, δτ=0.25, bands=2, contour=:imag, symmetry=sym)
            kernel = correlationfunction(bath, lat; band=1)
            IF1, _ = hybriddynamics(lat, kernel, PartialIF())
            vac = vacuumstate(lat)
            K = sysdynamics(lat, AndersonIM(U=U, ϵ_d=μ); order=2, coupling_norm=cn)
            zs = integrate_scaled(K, IF1, vac)
            Z = exp(zs.logabs) * zs.phase
            occ = occupation(lat, AndersonIM(U=U, ϵ_d=μ), nothing, IF1, vac; j=2, branch=:τ,
                             band=1, order=2, coupling_norm=cn, Z=Z)
            @test abs(Z - ref_Z) / abs(ref_Z) < 1e-10    # ③ Z 一致
            @test abs(occ - ref_occ) < 1e-10             # ③ 观测一致
        end
    end
    # ---------- 实时：演化 + 带因子化 ----------
    let U = 0.0, N = 3, dt = 0.15
        ed = ed_reference(U; kind=:real, dt=dt, Nt=N)
        lat = GrassmannLattice(N=N, δt=dt, bands=2, contour=:real)
        cn = α * 2
        kernel = correlationfunction(bath, lat)
        Is, _ = hybriddynamics(lat, kernel, PartialIF())
        ρ0 = Diagonal(ComplexF64[1.0, exp(-β * μ), exp(-β * μ), exp(-β * 2μ)]); ρ0 /= tr(ρ0)
        K = sysdynamics(lat, AndersonIM(U=0.0, ϵ_d=μ); order=2, coupling_norm=cn, initial_state=ρ0)
        zs = integrate_scaled(K, Is...)
        Z = exp(zs.logabs) * zs.phase
        # ① U=0 带因子化：⟨n₁⟩ ≡ M=1 值
        lat1 = GrassmannLattice(N=N, δt=dt, bands=1, contour=:real)
        IF1m, _ = hybriddynamics(lat1, correlationfunction(bath, lat1), PartialIF())
        K1 = sysdynamics(lat1, ToulouseIM(ϵ_d=μ); order=2, coupling_norm=α,
                                    initial_state=Diagonal(ComplexF64[1.0, exp(-β * μ)]) / (1 + exp(-β * μ)))
        zs1 = integrate_scaled(K1, IF1m)
        Z1 = exp(zs1.logabs) * zs1.phase
        occ_m2 = occupation(lat, AndersonIM(U=0.0, ϵ_d=μ), nothing, Is...; j=N, branch=:+, band=1,
                            order=2, coupling_norm=cn, initial_state=ρ0, Z=Z)
        occ_m1 = occupation(lat1, ToulouseIM(ϵ_d=μ), nothing, IF1m; j=N, branch=:+, band=1,
                            order=2, coupling_norm=α,
                            initial_state=Diagonal(ComplexF64[1.0, exp(-β * μ)]) / (1 + exp(-β * μ)), Z=Z1)
        @test abs(real(occ_m2) - real(occ_m1)) < 1e-10        # ① 因子化
        err = abs(real(occ_m2) - ed)
        @printf("  real U=0 N=3: ⟨n₁(t)⟩=%.7f  ED=%.7f  err=%.2e  (M=1≡M=2 ✓)\n", real(occ_m2), ed, err)
        @test err < 1e-4                                       # ② ED 锚
        # ④ U=0 带对称：双浴下 ⟨n₂⟩ = ⟨n₁⟩（同浴同模型逐位对称）
        occ2 = occupation(lat, AndersonIM(U=0.0, ϵ_d=μ), nothing, Is...; j=N, branch=:+, band=2,
                          order=2, coupling_norm=cn, initial_state=ρ0, Z=Z)
        @test abs(real(occ2) - real(occ_m2)) < 1e-10
    end
    println("Phase 5c 物理电池: all green")
end
