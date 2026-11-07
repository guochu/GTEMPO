# =============================================================================
# imag_pipeline_calib.jl —— GMPS 管线标定 v2（用稠密基准做判决锚）
#
# proto_imag.jl 的 24 候选全灭发生在 Iset bug 阶段（D 矩阵曾是 4×4 全槽位
# 错矩阵）；修复后的"纯 Δ = −D_full × siteflip × counter_c"从未测过。
# 本脚本：正确 D + 稠密真值（imag_dense.jl，机器精度）当判决锚，小网格
# 锁定 Δ⊕K 复合约定。
#
# 结构翻译（稠密 ↔ 管线，已验证稠密侧）：
#   稠密位置 p 的 (bar,eta) ↔ 管线切片 p+1 的 (ā,a)（corr 节点 p）；
#   稠密 state p 分辨率 (1−ξ̄_pξ_{p+1}) ↔ 管线切片 p+1 局部分辨率；
#   稠密区间 kernel(ξ̄_p,ξ_p) ↔ 管线门 g（切片 p→p+1）；
#   管线多一片 slice 1（冗余，BC exp(−ā₁a_k) 实现切片 k ≡ −切片 1 认同）。
#
# 内存：全小（N ≤ 3，稠密 L ≤ 16，GMPS NoTruncation 小键）——远低于闸。
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_pipeline_calib.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_dense.jl"))  # 稠密基准（include 安全）
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

# imag_dense 提供：BETA HB VB EPS_D hcoupled interval_operators dense_ring
#                 fock_ring ed_n；formb_dense 提供 Monomials 机器

# ------------------------------------------------- 环 M-Schur（正确 Iset）--
"""未拆分环 M 的 Schur 消元（proto_imag 修复版：杂质模 = slot s=1）。

返回 (Δ = −D_full, ν_a = det(M)/det(Mb), dMb)。"""
function ring_delta(cells::Vector{Float64})
    N = length(cells)
    h = hcoupled()
    mm = size(h, 1)
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0                 # 反周期缝（M[N,N] = +U_N）
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0  # 迹 resolution（前向对角）
        end
    end
    Iset = [(r-1)*mm + 1 for r in 1:N]              # 杂质模 = slot 1（修复点）
    Bset = [(r-1)*mm + s for r in 1:N for s in 2:mm]
    D = M[Iset, Iset] - M[Iset, Bset] * (M[Bset, Bset] \ M[Bset, Iset])
    Mb = zeros(ComplexF64, (mm-1)*N, (mm-1)*N)
    hbb = exp(-(BETA / N) .* HB)
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
    return -D, det(M) / det(Mb), det(Mb)
end

# ------------------------------------------------------- ImagDelta 适配器 --
mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
    siteflip::Bool
end
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
function Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.branch == :τ && j.branch == :τ) || throw(ArgumentError("branch must be :τ"))
    v = x.Δ[i.j, j.j]
    x.siteflip && i.j > j.j && (v = -v)   # partialmpo 位点序：行腿在列腿后翻号
    return v
end

# ------------------------------------------------------------- 管线运行 --
"""GMPS τ 环管线：环调度门 + IF(adapter Δ) + 接口反测度 cc + BC 折返。

返回 (Z_pipe, n)。Z_pipe = integrate·ν（ν 传入）。"""
function run_pipeline(Δ, ν, siteflip, cc, eps_d; N=2, order=2)
    sched = NativeExact.native_schedule(order, N, BETA / N)
    links = sched.links
    ring_links = vcat(links[1] + links[end], links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(Δ, siteflip); band=1, trunc=NoTruncation())
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    model = ToulouseIM(ϵ_d=eps_d)
    for g in 1:nl
        fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
        sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
        mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
    end
    if cc != 0
        for node in 2:(lat.k - 1)
            s_eta = lat[GrassmannContourIndex(node, :τ, 1, false)]
            s_bar = lat[GrassmannContourIndex(node, :τ, 1, true)]
            fac = GTEMPO.partialmpo(s_eta, [s_bar], [ComplexF64(cc)]) *
                  GTEMPO.GrassmannMPS(ComplexF64, length(lat))
            mpsK, _ = GTEMPO.mult(mpsK, fac; trunc=NoTruncation())
        end
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    Z = GTEMPO.integrate(lat, mpsK, mpsI) * ν
    b = GrassmannContourIndex(1, conj=true, branch=:τ, band=1)
    gv = GTEMPO.gf(lat, (GrassmannContourIndex(lat.Nτ, conj=false, branch=:τ, band=1), b),
                   mpsK, mpsI, Z=Z)
    return Z, real(gv)          # gτ(β⁻) 约定 = ⟨n⟩（firstorder 教程同款）
end

# ================================================================== 主程序 --
function main()
    for N in (2, 3)
        cells = fill(BETA / N, N); durs = copy(cells)
        # 稠密真值（判决锚）
        Ws, n1 = interval_operators(cells, durs, EPS_D)
        Wsn = [p == 1 ? n1 * Ws[1] : Ws[p] for p in 1:N]
        zd, _ = dense_ring(Ws)
        zdn, _ = dense_ring(Wsn)
        n_dense = real(zdn / zd)
        # 管线候选
        Δfull, ν_a, _ = ring_delta(cells)
        zref = real(zd)                                    # 物理链迹（未归一）
        println("─"^74)
        @printf("N=%d  稠密锚: Z = %.12f  n = %.12f\n", N, zref, n_dense)
        for sgn in (1, -1), siteflip in (false, true), cc in (0, 1, -1)
            Δ = sgn == 1 ? Δfull : -Δfull
            z, n = run_pipeline(Δ, ν_a, siteflip, cc, EPS_D; N=N)
            @printf("  sgn=%+d sf=%d cc=%+d : Z = %14.10f (err %9.2e)   n = %12.10f (err %9.2e)\n",
                    sgn, siteflip, cc, real(z), abs(real(z) - zref),
                    n, abs(n - n_dense))
        end
    end
end
main()
