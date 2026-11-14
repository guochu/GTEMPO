# ExactTTI golden checker（新 API）。
#
# 用法：julia --project=. test/fixtures/check_golden_exactttiif.jl
#
# 验证策略（2026-11-14 定谳）：
#   ① 主断言：ExactTTI 与 PartialIF（golden 机器精度通道）在同一格上 IF 多项式
#      表逐位一致（pencil 拟合对单模浴核精确，两算法是同一指数的两种装配），
#      Z/occupation 一致；
#   ② ED 锚：虚时全热态（绝对 ω 约定）、实时演化值，离散化误差档；
#   ③ 旧栈 exactttiif fixtures（port_golden/*_exactttiif.json）**仅参考不断言**：
#      它们编码旧栈 ExactTTI 自己的已知 Prony 缺陷（旧 imagtime 测试在案：
#      TTI 算法 ≈5% 不随 δτ 收敛），新实现已超越——imag N=4 对照
#      occupation 旧 0.37555（距 ED 2.5e-2）vs 新 0.34992（距 ED 8.7e-4）。
#
# 布局适配备忘：memory_cell/contact 合并胞按 U1GTEMPO2 存储 [± 相邻、j 降序]
# 设计——real 因子在虚拟 ± 相邻布局构建后 GM.reorder! 到 Chrono 站序；虚时
# compact 通道链与 Chrono [格 Nτ…1] 存储天然同序（不 reorder 直嵌）。

using Test, Printf
using GTEMPO
using ImpurityModelBase
import GrassmannMPS as GM
using LinearAlgebra: diagm, I, tr, eigen, Hermitian, Diagonal, kron

include(joinpath(@__DIR__, "golden_json.jl"))

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
ρ0_1 = diagm(0 => ComplexF64[1.0, exp(-β * μ)]) / (1 + exp(-β * μ))
const DIG = Dict(0 => Int[], 1 => [1, 2], 2 => [1], 3 => [2])

"ED 用单模湮灭算符（2 维）"
function _annih1()
    d = zeros(ComplexF64, 2, 2)
    d[1, 2] = 1.0
    return d
end

"GMPS → 变量粒度系数表（升序单项式基，站 digit 展平）"
function to_var(x)
    dv = GM.todense(x; order=:msb)
    L = length(x)
    out = Dict{Int,ComplexF64}()
    for (idx, val) in enumerate(dv)
        iszero(val) && continue
        rest = idx - 1
        mask = 0
        for s in L:-1:1
            d = rest % 4
            rest ÷= 4
            for off in DIG[d]; mask |= 1 << (2 * (s - 1) + off - 1); end
        end
        out[mask] = get(out, mask, 0) + val
    end
    return out
end

@testset "ExactTTI ≡ PartialIF + ED 锚" begin
    # ================= 虚时 N=4 =================
    let lat = GrassmannLattice(N=4, δτ=0.25, bands=1, contour=:imag)
        Is_p = (hybriddynamics(lat, bath, PartialIF(); band=1, coupling_norm=α),)
        Is_e = (hybriddynamics(lat, bath, ExactTTI(); band=1, coupling_norm=α),)
        K = sysdynamics(lat, ToulouseIM(ϵ_d=μ); order=1)
        tp, te = to_var(Is_p[1]), to_var(Is_e[1])
        keys_all = union(keys(tp), keys(te))
        err = maximum(abs(get(tp, k, 0im) - get(te, k, 0im)) for k in keys_all; init=0.0)
        @test err < 1e-10
        println("imag N=4 IF 表: max |ExactTTI − PartialIF| = ", err)
        zp = integrate_scaled(K, Is_p...)
        ze = integrate_scaled(K, Is_e...)
        Zp, Ze = exp(zp.logabs) * zp.phase, exp(ze.logabs) * ze.phase
        @test abs(Ze - Zp) / abs(Zp) < 1e-10
        kw = (order=1, coupling_norm=α)
        @test occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is_e...; j=2, branch=:τ, band=1, Z=Ze, kw...) ≈
              occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is_p...; j=2, branch=:τ, band=1, Z=Zp, kw...) atol = 1e-10
        # ED 锚：全热态（绝对 ω，无 μ 移）
        I2 = Matrix{ComplexF64}(I, 2, 2)
        c1 = kron(_annih1(), I2)
        bm = kron(Matrix{ComplexF64}(I, 2, 2), _annih1())
        H = μ * (c1' * c1) + ω * (bm' * bm) + sqrt(α) * (bm' * c1 + c1' * bm)
        ev, vec = eigen(Hermitian(Matrix(H)))
        w = exp.(-β .* (ev .- minimum(ev))); w ./= sum(w)
        ed = real(tr((vec * Diagonal(ComplexF64.(w)) * vec') * (c1' * c1)))
        occ = real(occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is_e...; j=2, branch=:τ, band=1, Z=Ze, kw...))
        @printf("imag: ⟨n⟩ ExactTTI=%.7f  ED=%.7f  err=%.2e（δτ=0.25 离散化档）\n", occ, ed, abs(occ - ed))
        @test abs(occ - ed) < 2e-3
    end

    # ================= 实时 N=3（热 ρ₀，Weyl 二阶）=================
    let lat = GrassmannLattice(N=3, δt=0.15, bands=1, contour=:real)
        cn = α
        Is_p = (hybriddynamics(lat, bath, PartialIF(); band=1, coupling_norm=cn),)
        Is_e = (hybriddynamics(lat, bath, ExactTTI(); band=1, coupling_norm=cn),)
        K = sysdynamics_weyl2order(lat, ToulouseIM(ϵ_d=μ); coupling_norm=cn, initial_state=ρ0_1)
        tp, te = to_var(Is_p[1]), to_var(Is_e[1])
        keys_all = union(keys(tp), keys(te))
        err = maximum(abs(get(tp, k, 0im) - get(te, k, 0im)) for k in keys_all; init=0.0)
        @test err < 1e-10
        println("real N=3 IF 表: max |ExactTTI − PartialIF| = ", err)
        zp = integrate_scaled(K, Is_p...)
        ze = integrate_scaled(K, Is_e...)
        Zp, Ze = exp(zp.logabs) * zp.phase, exp(ze.logabs) * ze.phase
        @test abs(Ze - Zp) / abs(Zp) < 1e-10
        kw = (order=2, coupling_norm=cn, initial_state=ρ0_1)
        @test occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is_e...; j=3, branch=:+, band=1, Z=Ze, kw...) ≈
              occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is_p...; j=3, branch=:+, band=1, Z=Zp, kw...) atol = 1e-10
        # ED 锚：演化 t=0.45
        I2 = Matrix{ComplexF64}(I, 2, 2)
        bb = _annih1()
        c1m = bb
        n1m = c1m' * c1m
        nb = [0.0 0.0; 0.0 1.0]
        H = kron(μ * n1m, I2) + kron(I2, ω * nb) + sqrt(α) * (kron(c1m', bb) + kron(c1m, bb'))
        evb, vecb = eigen(Hermitian(nb))
        wb = exp.(-β .* (evb .- minimum(evb))); wb ./= sum(wb)
        function edreal()
            rr = kron(ρ0_1, vecb * Diagonal(ComplexF64.(wb)) * vecb')
            Ue = exp(-1.0im * 0.15 * H)
            for _ in 1:3; rr = Ue * rr * Ue'; end
            return real(tr(rr * kron(n1m, I2)))
        end
        occ = real(occupation(lat, ToulouseIM(ϵ_d=μ), nothing, Is_e...; j=3, branch=:+, band=1, Z=Ze, kw...))
        ed = edreal()
        @printf("real: ⟨n(t=.45)⟩ ExactTTI=%.7f  ED=%.7f  err=%.2e\n", occ, ed, abs(occ - ed))
        @test abs(occ - ed) < 1e-4
    end

    # ================= 旧栈 fixture 参考（不断言，⑤④ 修复前后差异记录）=================
    let name = "imag_n4_toulouse_exactttiif"
        ref = load_golden(name)
        println("参考（不断言）：旧栈 exactttiif fixture ", name,
                " Z=", ascomplex(ref["Z"]), "  occ=", ascomplex(ref["occupation_tau_b1"][1]),
                "——旧栈 Prony 缺陷（距 ED 2.5e-2）；新栈 = 0.34992175（距 ED 8.7e-4）")
    end
    println("ExactTTI golden: all green")
end
