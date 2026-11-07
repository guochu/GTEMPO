# =============================================================================
# imag_calib_fit.jl —— Δ⊕K 复合约定锁定：真机器 oracle + Gauss–Newton 拟合
#
# 背景（imag_skeleton_extract.jl / imag_trace_probe.jl 的结论）：
#   · 提取约定已锁定（vacuum/BC/IF-分析式 全机器精度）；
#   · Δ 的 pair 直读 = D_full（−Δfull），但端到端 Z/n 未过；
#   · integrate 的 trace 配对是"每切片 K/I 二选一分辨率"型奇特泛函
#     （17 探针自洽），跨界线对取值未定——反推代数性价比低。
# 本脚本：把真实 run_pipeline 当 oracle，拟合 Δ（实参数化，h 实 ⇒ D 实），
# 残差 = Z(εd 网格) vs Fock 环锚。初始点 D_full。泛化判决 = order 1 / N=3。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_calib_fit.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

# ------------------------------------------------------------------ 模型 --
const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

# ------------------------------------------------- 环 M-Schur（正确 Iset）--
function ring_delta(cells::Vector{Float64})
    N = length(cells)
    h = hcoupled()
    mm = size(h, 1)
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0
        end
    end
    Iset = [(r-1)*mm + 1 for r in 1:N]
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
    return -D, det(M) / det(Mb)
end

# ------------------------------------------------------- ImagDelta 适配器 --
mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
end
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex) = x.Δ[i.j, j.j]

# ------------------------------------------------------------- 管线运行 --
"""GMPS τ 环管线。返回 (Z = integrate·ν, gτ 变体组)。"""
function run_pipeline(Δ, ν, eps_d; N=2, order=2)
    sched = NativeExact.native_schedule(order, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(Δ); band=1, trunc=NoTruncation())
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    model = ToulouseIM(ϵ_d=eps_d)
    for g in 1:nl
        fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
        sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
        mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    Z = GTEMPO.integrate(lat, mpsK, mpsI) * ν
    b1 = GrassmannContourIndex(1, conj=true, branch=:τ, band=1)
    g_end = real(GTEMPO.gf(lat, (GrassmannContourIndex(lat.Nτ, conj=false, branch=:τ, band=1), b1),
                           mpsK, mpsI, Z=Z))
    g_1 = real(GTEMPO.gf(lat, (GrassmannContourIndex(1, conj=false, branch=:τ, band=1), b1),
                         mpsK, mpsI, Z=Z))
    return Z, (g = g_end, g1 = g_1, om1 = 1 - g_1, ome = 1 - g_end)
end

"""Fock 环锚：Z 与 n（dim = 2·2^m 多体，门在 cell 后作用链）。"""
annih_op(dim::Int, mode::Int) = begin          # formb_dense 同款（JW 符号）
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    return op
end
function fock_ring(cells, gate_durs, eps_d)
    m = length(VB)
    nt = 1 + m
    dim = 1 << nt
    ops = [annih_op(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) && continue
        Hq .+= h[i, j] .* (ops[i]' * ops[j])
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Z = Matrix{ComplexF64}(I, dim, dim)
    Zn = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z
        Zn = W * Zn
    end
    Zn = ops[1]' * ops[1] * Zn
    return real(tr(Z)), real(tr(Zn) / tr(Z))
end

# ---------------------------------------------------------- Gauss–Newton --
"""实参数化 Δ(x)：x ∈ R^{N²} → Δ 矩阵。"""
function delta_of(x::Vector{Float64}, N::Int)
    Δ = zeros(ComplexF64, N, N)
    for j in 1:N, i in 1:N
        Δ[i, j] = x[(j-1)*N + i]
    end
    return Δ
end

"""残差向量：Z(εd_j) − Z锚_j（order 2, N=2, ν=1）。"""
function residuals(x, grid, zanchors; N=2, order=2, ν=1.0)
    Δ = delta_of(x, N)
    return [real(run_pipeline(Δ, ν, e; N=N, order=order)[1]) - z for (e, z) in zip(grid, zanchors)]
end

function gauss_newton(x0, grid, zanchors; N=2, order=2, ν=1.0, iters=12, damp=1e-6)
    x = copy(x0)
    for it in 1:iters
        r = residuals(x, grid, zanchors; N=N, order=order, ν=ν)
        J = zeros(length(r), length(x))
        for k in eachindex(x)
            xp = copy(x); xp[k] += 1e-7
            xm = copy(x); xm[k] -= 1e-7
            J[:, k] = (residuals(xp, grid, zanchors; N=N, order=order, ν=ν) -
                       residuals(xm, grid, zanchors; N=N, order=order, ν=ν)) ./ 2e-7
        end
        dx = (J' * J + damp * I) \ (-J' * r)
        x += dx
        @printf("  iter %2d : |r| = %9.3e  |dx| = %9.3e\n", it, norm(r), norm(dx))
        norm(r) < 1e-12 && break
        norm(dx) < 1e-14 && break
    end
    return x, norm(residuals(x, grid, zanchors; N=N, order=order, ν=ν))
end

# ================================================================== 主程序 --
function main()
    N = 2
    cells = fill(BETA / N, N)
    sched = NativeExact.native_schedule(2, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    Δfull, ν_a = ring_delta(cells)

    grid = [0.0, 0.35, 0.7, 1.05, 1.4]
    zanchors = [fock_ring(cells, ring_links, e)[1] for e in grid]
    nanchors = [fock_ring(cells, ring_links, e)[2] for e in grid]
    println("锚（order 2, N=$N, cells=$cells, links=$ring_links）：")
    for (e, z, n) in zip(grid, zanchors, nanchors)
        @printf("  εd=%.2f : Z = %.12f  n = %.12f\n", e, z, n)
    end
    @printf("D_full = "); display(round.(Δfull, digits=6))
    @printf("ν_a = %.12f\n", real(ν_a))

    println("\n── D_full 初值的残差 ──")
    x0 = vec(real.(Δfull))
    r0 = residuals(x0, grid, zanchors; N=N)
    @printf("x0 = %s  |r| = %.3e\n", [round(v, digits=4) for v in x0], norm(r0))

    println("\n── Gauss–Newton 拟合（实 4 参数，ν=1）──")
    xfit, rf = gauss_newton(x0, grid, zanchors; N=N)
    @printf("拟合 Δ = \n"); display(round.(delta_of(xfit, N), digits=8))
    @printf("|r_final| = %.3e\n", rf)

    println("\n── 拟合 Δ 的 n 判决（gτ 变体 vs 锚）──")
    Δfit = delta_of(xfit, N)
    for (e, na) in zip(grid[2:4], nanchors[2:4])
        _, g = run_pipeline(Δfit, 1.0, e; N=N)
        @printf("εd=%.2f : 锚 n=%.8f | gτ(end)=%9.6f gτ(1)=%9.6f 1−gτ(1)=%9.6f 1−gτ(end)=%9.6f\n",
                e, na, g.g, g.g1, g.om1, g.ome)
    end

    println("\n── 泛化判决 1：order 1（同一 Δ_fit，不同 cells/links）──")
    sched1 = NativeExact.native_schedule(1, N, BETA / N)
    rl1 = vcat(sched1.links[1] + sched1.links[end], sched1.links[2:end-1])
    Δf1, _ = ring_delta(fill(BETA / N, N))   # order 1: cells 同
    z1 = fock_ring(fill(BETA / N, N), rl1, 0.7)[1]
    zp1, g1 = run_pipeline(Δfit, 1.0, 0.7; N=N, order=1)
    @printf("order 1 εd=0.7 : Z 管线 %.10f vs 锚 %.10f  err %.2e\n", real(zp1), z1, abs(real(zp1) - z1))

    println("\n── 泛化判决 2：N=3（重新 ring_delta，Δ_fit 仅作起点参考）──")
    N3 = 3
    cells3 = fill(BETA / N3, N3)
    sched3 = NativeExact.native_schedule(2, N3, BETA / N3)
    rl3 = vcat(sched3.links[1] + sched3.links[end], sched3.links[2:end-1])
    Δfull3, _ = ring_delta(cells3)
    z3 = fock_ring(cells3, rl3, 0.7)[1]
    zp3, g3 = run_pipeline(Δfull3, 1.0, 0.7; N=N3)
    @printf("N=3 D_full  εd=0.7 : Z = %.10f vs 锚 %.10f  err %.2e\n", real(zp3), z3, abs(real(zp3) - z3))
    _, g3d = run_pipeline(Δfull3, 1.0, 0.7; N=N3)
    @printf("N=3 锚 n=%.8f | gτ(end)=%.6f 1−gτ(end)=%.6f\n", fock_ring(cells3, rl3, 0.7)[2], g3d.g, g3d.ome)
end
main()
