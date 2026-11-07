# =============================================================================
# proto_imag.jl —— 严格 IF 虚轴原型：环 Gaussian 消元（已验证）+ 管线标定（进行中）
#
# 【已锁定（2026-11-07）】
# 1. 环 Gaussian 装配精确：相干态迹网格公式
#      ⟨ξ̄_p|O_p|ξ_p⟩ = exp(ξ̄_pU_pξ_p) → M[p,p] = −U_p（p<N），缝位反周期
#      ⟨-ξ̄_N| ⇒ M[N,N] = +U_N；resolution 𝟙 = ∫|ξ⟩⟨ξ̄|e^{-ξ̄ξ}（每插入一个）
#      → M[p,p+1] = +1（p<N）、缝 resolution M[N,1] = +1。
#    det(M) == Tr_多体[e^{-βh₀ 链}] == det(I+∏U)（Krein）机器精度（N=2,3,4 验证；
#    符号因子 (−1)^{(1+m)(N-1)}，bath-only 同构环同因子 ⇒ ν = det(M)/det(Mb) 无符号）。
#    （曾犯的错：resolution 写成同位对角 M[p,p] += 1——错；`1 .+ Matrix` 逐元素加
#    ≠ I+A——Krein 检验必须 det(I+U)；Iset 须只含杂质模 slot s=1。）
# 2. Schur 消 bath：D_full = M_II − M_IB M_BB⁻¹ M_BI（含 resolution）；D_cell =
#    去杂质模 resolution 对角后的同式。归一化候选 ν_a = det(M)/det(Mb)。
# 3. 骨架（环调度门：首末半门环合并 → nl 门；GTEMPO boundarycondition! 折返）
#    在 Δ=0 时精确给出自由杂质环迹 Z₀ = 2（εd=0）——骨架测度干净、无多余因子。
# 4. 骨架在 εd=0 是**高斯**多项式（Z(Δ) 探测线性响应干净）：
#    探测 P[i,j] = ∂Z/∂Δᵢⱼ / Z₀ 呈闭式结构 P[i,j] = +0.5 (i>j)、−0.5 (i≤j)。
#    注意 P 是 f(K₁)（= det(K₁+Δ)/det(K₁) 型配分函数）的对数导数，不是骨架
#    双线性矩阵 K₁ 本身。
#
# 【未解的结（下一步）】
# 管线隐式测度配对（integrate/gf 对 imag lattice 的 Grassmann 收缩语义，
# 含 slice 1 非_IF 腿与 BC 折返项）与环槽位（corr 节点 i ↔ lattice 切片 i+1）
# 的复合约定未锁定：候选 Δ = −D_full/D_cell − K₁（P 的 invP / invPT 解释，
# 各 × siteflip 开关）共 24 组合无一命中 n(β)==ED。需要 formb_dense Step B/C
# 式的**稠密多项式副本**：把骨架（门张量 + BC + 隐式 resolution）展开为显式
# Grassmann 多项式、与环消元的杂质扇区泛函做逐系数 diff，直接读出 Δ 与 K
# 的复合映射（含 partialmpo 位点序符号——实时路线 siteflip 的虚轴对应物）。
#
# 本文件现状：imag_ring_delta（装配+消元，已验证）、ImagDelta 适配器、
# run_pipeline（环管线 + n(β) 读出）、probe_main（K₁ 探测 + 候选检验）、
# grid_main（旋钮网格，已被证伪保留作记录）。入口 = probe_main()。
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation

include(joinpath(@__DIR__, "..", "src", "NativeExact.jl"))
using .NativeExact

# ------------------------------------------------------------------ 模型 --
const BETAP = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]          # DiracDelta ⇔ |V|² = α

"""耦合单粒子 h（复用 _coupled_h：mode 1 = 杂质 ε=0）。"""
hcoupled() = NativeExact._coupled_h(VB, HB)

"""环 Gaussian 消元 → (Δ, ν)。

相干态迹网格（已验证 N=1: det(1+U)；N=2: (−1)^{N-1} det(1+U₂U₁)）：
- ⟨ξ̄_p|O_p|ξ_p⟩ = exp(ξ̄_p U_p ξ_p)：M[p,p] = −U_p（p<N），缝位反周期
  ⟨-ξ̄_N| ⇒ M[N,N] = +U_N；
- resolution 𝟙 = ∫|ξ⟩⟨ξ̄|e^{-ξ̄ξ}（每插入一个）：M[p, p+1] = +1（p<N），
  缝 resolution M[N,1] = +1。
det(M) = (−1)^{(1+m)(N-1)}·det(1 + U_N···U_1)；ν = det(M)/det(Mb)（同构 bath
环，符号因子相消）= det(1+e^{-βh₀ 链})/det(1+e^{-βh_bath})。"""
function imag_ring_delta(cells::Vector{Float64})
    N = length(cells)
    h = hcoupled()
    mm = size(h, 1)               # 1 + m
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0        # 反周期缝
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0    # resolution：bra_r → ket_{r+1}
        end
    end
    ii = 1                        # 杂质模 = slot s=1（bath = 2..mm）
    Iset = [ (r-1)*mm + s for r in 1:N for s in ii ]
    Bset = [ (r-1)*mm + s for r in 1:N for s in 2:mm ]
    MII = M[Iset, Iset]; MIB = M[Iset, Bset]
    MBB = M[Bset, Bset]; MBI = M[Bset, Iset]
    D = MII - MIB * (MBB \ MBI)
    # 候选 D_cell：不含杂质模 resolution 对角的环（resolution 由骨架 K 承载；
    # bath 侧 resolution 属于本消元积分，保留）
    Mc = copy(M)
    for r in 1:N
        nxt = r == N ? 1 : r + 1
        Mc[(r-1)*mm + 1, (nxt-1)*mm + 1] -= 1.0
    end
    Ic = Mc[Iset, Iset]; IBc = Mc[Iset, Bset]; BIc = Mc[Bset, Iset]
    Dcell = Ic - IBc * (MBB \ BIc)
    # bath-only 同构环（归一化参照；符号因子与 det(M) 相消）
    Mb = zeros(ComplexF64, (mm-1)*N, (mm-1)*N)
    hbb = exp(-(BETAP / N) .* HB)    # bath-only cell = e^{-cell·h_bath}（均匀 cells）
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm-1, s1 in 1:mm-1
            Mb[(r-1)*(mm-1) + s2, (r-1)*(mm-1) + s1] += sgn * hbb[s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm-1
            Mb[(r-1)*(mm-1) + s, (nxt-1)*(mm-1) + s] += 1.0
        end
    end
    return -D, Dcell, det(M), det(Mb), det(MBB), Mb
end

"""单粒子独立性检验：det(M) vs (−1)^k·det(1 + U_N···U_1)。"""
function check_ring_det(cells)
    h = hcoupled()
    N = length(cells)
    U = Matrix{ComplexF64}(I, size(h,1), size(h,1))
    for r in N:-1:1
        U = exp(-cells[r] .* h) * U
    end
    return det(I + U)
end

# ------------------------------------------------------- ImagDelta 适配器 --
mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
    ν::ComplexF64
    sgn::Int
    siteflip::Bool
end
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
function Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.branch == :τ && j.branch == :τ) || throw(ArgumentError("branch must be :τ"))
    (i.band == j.band == 1) || throw(ArgumentError("single band only in proto"))
    v = x.Δ[i.j, j.j]
    x.sgn == 1 || (v = -v)
    # partialmpo 位点序：行腿（bar，切片 i+1）在列腿（eta，切片 j+1）之后时翻号
    x.siteflip && i.j > j.j && (v = -v)
    return v
end
branch(x::ImagDelta, b1::Symbol, b2::Symbol) = ((b1 == :τ) && (b2 == :τ)) || throw(ArgumentError("branch must be :τ"))

# ------------------------------------------------------------ 管线运行 --
function run_pipeline(Δ, ν, sgn, counter_c, eps_d; N=2, order=1, want_n=true, siteflip=false)
    sched = NativeExact.native_schedule(order, N, BETAP / N)
    links = copy(sched.links)
    # 环合并：首末半门并成一个门（放在 interval 1）
    ring_links = vcat(links[1] + links[end], links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETAP / nl, contour=:imag, bands=1)
    d = ImagDelta(Δ, ν, sgn, siteflip)
    d.ν = ν
    mpsI = GTEMPO.hybriddynamics(lat, d; band=1, trunc=NoTruncation())
    # 骨架：环调度门 + 反测度（切片 2..k）+ BC
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    model = ToulouseIM(ϵ_d=eps_d)
    for g in 1:nl
        fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
        sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
        mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
    end
    for node in 2:(lat.k - 1)     # 接口反测度候选（缝位 1 与 k 不放）
        s_eta = lat[GrassmannContourIndex(node, :τ, 1, false)]
        s_bar = lat[GrassmannContourIndex(node, :τ, 1, true)]
        fac = GTEMPO.partialmpo(s_eta, [s_bar], [ComplexF64(counter_c)]) *
              GTEMPO.GrassmannMPS(ComplexF64, length(lat))
        mpsK, _ = GTEMPO.mult(mpsK, fac; trunc=NoTruncation())
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    Z = GTEMPO.integrate(lat, mpsK, mpsI) * ν
    want_n || return Z, 0.0
    # ⟨n(β)⟩：缝位折返腿 (a⁺_k, ā⁻_k) 上插 n̂ = a†a → GTerm(ā_k, a_k) 插入
    b = GrassmannContourIndex(1, conj=true, branch=:τ, band=1)   # corr 节点 1 = 切片 2
    nvals = [GTEMPO.gf(lat, (GrassmannContourIndex(i, conj=false, branch=:τ, band=1), b),
                       mpsK, mpsI, Z=Z) for i in 1:lat.Nτ]
    # gτ(i) = ⟨a(τᵢ)a†(τ₁)⟩；n(β) = 1 − gτ(N)（τ=β⁻ 反周期像）
    return Z, real(1 - nvals[end])
end

# ---------------------------------------------------------------- ED 参照 --
function ed_thermal(eps_d, beta)
    dim = 4
    a1 = zeros(ComplexF64, dim, dim); a2 = zeros(ComplexF64, dim, dim)
    a1[2, 1] = 1.0; a2[4, 3] = 1.0
    n1 = a1' * a1; n2 = a2' * a2
    H = conj(VB[1]) * (a1' * a2) + VB[1] * (a2' * a1) .+ HB[1,1] .* n2 .+ eps_d .* n1
    R = exp(-beta .* H)
    return R / tr(R), n1, H
end

# ================================================================== 主程序 --

# ------------------------------------------------ K₁ 探测（骨架双线性结构）--
"""εd=0 时 K = Z₀·exp(ξ̄K₁ξ)（门=恒等 ⇒ 骨架纯测度、高斯）。
对 Δ = t·Eᵢⱼ 有限差分提 K₁[i,j]。（ν=1, sgn=+1, cc=0 档）"""
function probe_K1(N::Int)
    cells = fill(BETAP / N, N)
    Δfull, Dc, dM, dMb, dMBB, _ = imag_ring_delta(cells)
    Z0, _ = run_pipeline(zero(Δfull), 1.0 + 0.0im, 1, 0, 0.0; N=N, want_n=false)
    t = 1.0e-4
    K1 = zeros(ComplexF64, N, N)
    for i in 1:N, j in 1:N
        E = zeros(ComplexF64, N, N); E[i, j] = 1.0
        zp, _ = run_pipeline(t * E, 1.0 + 0.0im, 1, 0, 0.0; N=N, want_n=false)
        zm, _ = run_pipeline(-t * E, 1.0 + 0.0im, 1, 0, 0.0; N=N, want_n=false)
        K1[i, j] = (zp - zm) / (2t * Z0)
    end
    return K1, Z0, Δfull, Dc, dM, dMb
end

function probe_main()
    for N in (2, 3, 4)
        K1, Z0, Δfull, Dc, dM, dMb = probe_K1(N)
        println("N=$N  Z(Δ=0)=Z₀ = ", real(Z0))
        println("K₁ ="); display(round.(K1, digits=6))
        # 候选：Δ_final = −D − K₁actual（探测给 f(K₁) 的对数导数 = K₁⁻¹ 转置项）
        K1inv = inv(K1)
        for (nm, Dbasis) in (("full", Δfull), ("cell", Dc)), (kn, Kact) in (("invP", K1inv), ("invPT", transpose(K1inv)))
            _, n1op, Hfull = ed_thermal(0.7, BETAP)
            ned = real(tr(n1op * exp(-BETAP .* Hfull)) / tr(exp(-BETAP .* Hfull)))
            for sf in (false, true)
                Z, n = run_pipeline(-Dbasis - Kact, dM / dMb, 1, 0, 0.7; N=N, siteflip=sf)
                @printf("  Δ=-D_%s-K₁(%s) sf=%d: n = %9.6f（ED %.6f, err %.2e）  Z = %10.6f\n",
                        nm, kn, sf, n, ned, abs(n - ned), real(Z))
            end
        end
    end
end

function grid_main()
    # 0) det 恒等式检验
    for N in (2, 3, 4)
        cells = fill(BETAP / N, N)
        Δ, Dc, dM, dMb, dMBB, Mb = imag_ring_delta(cells)
        zsp = check_ring_det(cells)     # det(I + ∏U) 单粒子
        @printf("N=%d  det(M)/%s = %.12f  det(1+∏U) = %.12e  det(MBB) = %.6e  det(Mb) = %.6e\n",
                N, "(−1)^k", dM, zsp, dMBB, dMb)
        @printf("      符号因子检验: det(Mb)/(bath det) = %.6f；det(M)/det(I+∏U) = %.6f\n",
                dMb / det(I + exp(-BETAP .* HB)), dM / zsp)
    end
    # 1) 物理锚点标定：n(β) vs ED（εd=0.7）
    _, n1op, Hfull = ed_thermal(0.7, BETAP)
    ned = real(tr(n1op * exp(-BETAP .* Hfull)) / tr(exp(-BETAP .* Hfull)))
    println("\nED: n(β) = ", ned)
    for N in (2, 3)
        cells = fill(BETAP / N, N)
        Δfull, Dc, dM, dMb, dMBB, _ = imag_ring_delta(cells)
        for (dname, Dbasis) in (("full", Δfull), ("cell", Dc)),
            sgn in (1, -1),
            (nuname, nu) in (("a", dM / dMb), ("b", dMBB / dMb), ("inv", dMb / dM)),
            cc in (0, 1, -1)
            Δ = sgn == 1 ? Dbasis : -Dbasis
            Z, n = run_pipeline(Δ, nu, 1, cc, 0.7; N=N)
            @printf("N=%d Δ=%-4s sgn=%+d ν=%-3s cc=%+d : n = %9.6f（err %9.2e）  Z = %10.5f\n",
                    N, dname, sgn, nuname, cc, n, abs(n - ned), real(Z))
        end
    end
end
probe_main()
# grid_main()
