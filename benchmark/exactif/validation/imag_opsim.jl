# =============================================================================
# imag_opsim.jl —— 算符层模拟器 + 注入链复合的标定探针（nl=2 小格穷举）
#
# 模型：imag-lattice 的 GMPS ↔ 切片副本 Fock 空间 H = ⊗_{s=1..k} C² 上的算符
#   （每切片 bar/eta = 该副本模式的产生/湮灭）。级联 apply!(exp(GTerm)) ↔
#   算符有序积；integrate(K, I) ↔ 某个双线性泛函 𝒯(O_K, O_I)。
# 标定目标：用 GTEMPO 实测复合量，扫描约定组合（嵌入方向 × 级联序 × 泛函）
#   找到全探针匹配的模型。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_opsim.jl [nl]
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0

# ------------------------------------------------------------ 模式算符 --
"""k 模式 Fock 空间（dim 2^k）上的 c_s / c†_s（JW 串跨模式）。基 |n_1…n_k⟩ 位掩码。"""
function mode_ops(k::Int)
    dim = 1 << k
    cs = [zeros(ComplexF64, dim, dim) for _ in 1:k]
    for n in 0:(dim-1)
        for s in 1:k
            ((n >> (s-1)) & 1) == 1 || continue
            jw = isodd(count_ones(n & ((1 << (s-1)) - 1))) ? -1.0 : 1.0
            m = n ⊻ (1 << (s-1))
            cs[s][m + 1, n + 1] = jw          # c|n⟩ = JW·|n⊖s⟩
        end
    end
    cds = [adjoint(cs[s]) for s in 1:k]
    return cs, cds
end

# ------------------------------------------------------------------ 主程序 --
function main()
    nl = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)                       # 2(nl+1)
    k = nl + 1
    @printf("nl=%d L=%d k=%d dim(H)=%d\n", nl, L, k, 1 << k)

    # 位图：slice s（1..k）× (bar,eta) → 1-based GTerm 位置
    pos = Dict{Tuple{Int,Bool}, Int}()
    for s in 1:k
        pos[(s, true)] = LAT[GrassmannContourIndex(s, :τ, 1, true)]    # bar
        pos[(s, false)] = LAT[GrassmannContourIndex(s, :τ, 1, false)]  # eta
    end
    println("位置表（1-based GTerm / 0-based pos）：")
    for s in 1:k
        @printf("  slice %d: bar %d/%d  eta %d/%d\n", s, pos[(s,true)], pos[(s,true)]-1,
                pos[(s,false)], pos[(s,false)]-1)
    end
    spos = Dict((s, b) => pos[(s, b)] - 1 for s in 1:k for b in (true, false))

    cs, cds = mode_ops(k)
    I_H = Matrix{ComplexF64}(I, 1 << k, 1 << k)
    vacuum0 = zeros(ComplexF64, 1 << k); vacuum0[1] = 1.0   # |0…0⟩

    """因子像：exp(GTerm(A 0-based, c)) ↔ I + c·∏_{p∈sort(A)} X_p。
    conv: bar→cd, eta→c（嵌入方向旋钮 barasc=false 时反向）。"""
    function factor_img(A0::Vector{Int}, c::ComplexF64; bar_to_creation=true)
        X = I_H * 0.0
        isfirst = true
        for p in sort(A0)
            s = findfirst(x -> x == p, [spos[(s, b)] for s in 1:k for b in (true, false)])
            # p → (slice, bar?) 查表
            sl = nothing; bar = nothing
            for ss in 1:k, bb in (true, false)
                spos[(ss, bb)] == p && (sl = ss; bar = bb)
            end
            sl === nothing && error("位置 $p 不属于任何切片腿")
            Xp = bar ? (bar_to_creation ? cds[sl] : cs[sl]) :
                       (bar_to_creation ? cs[sl] : cds[sl])
            X = isfirst ? Xp : X * Xp
            isfirst = false
        end
        return I_H + c .* X
    end

    """级联像：因子列表按 apply! 顺序相乘。order=:left（后者乘左）/:right。"""
    function cascade_img(factors::Vector{Tuple{Vector{Int},ComplexF64}}; order=:left)
        O = I_H
        for (A, c) in factors
            F = factor_img(A, c)
            O = order == :left ? F * O : O * F
        end
        return O
    end

    # ------------------------------------------------------- GTEMPO 侧 --
    EVENS = Vector{Vector{Int}}([Int[]])
    for m in 2:2:L, comb in combinations(0:(L-1), m)
        push!(EVENS, collect(comb))
    end
    injected = function(pairs::Vector{Tuple{Int,Int}})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (p1, p2) in sort(pairs)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mps)
        end
        return mps
    end
    inj_mps = Dict(A => injected([(A[2i-1], A[2i]) for i in 1:(length(A)÷2)]) for A in EVENS)
    vac = GTEMPO.vacuumstate(ComplexF64, LAT)

    meas(K, Ip) = real(GTEMPO.integrate(LAT, K, Ip))

    # --------------------------------------------------------- 探针电池 --
    println("\n== P0/P1/P2：真空侧行（τ2 表的 K/I 行）==")
    z00 = meas(vac, vac)
    println("P0 integrate(vac, vac) = ", z00)
    row1 = [meas(vac, inj_mps[A]) for A in EVENS]
    row2 = [meas(inj_mps[A], vac) for A in EVENS]
    println("P1 [integrate(vac, inj A)]（前 12）：", row1[1:min(12, end)])
    println("P2 [integrate(inj A, vac)]（前 12）：", row2[1:min(12, end)])

    # 模型候选扫描：嵌入方向 × 级联序（单因子平凡）× 泛函
    function model_val(OK::Matrix{ComplexF64}, OI::Matrix{ComplexF64}, fun::Symbol)
        fun === :vev_KI && return real((vacuum0' * (OK * OI) * vacuum0)[1])
        fun === :vev_IK && return real((vacuum0' * (OI * OK) * vacuum0)[1])
        fun === :tr && return real(tr(OK * OI))
        fun === :vev_Tr && return real((vacuum0' * OK * OI * vacuum0)[1] + (vacuum0' * OI * OK * vacuum0)[1])
        error("unknown fun")
    end
    println("\n== 约定扫描（P0 目标 = $z00）==")
    for bar2c in (true, false), ord in (:left, :right), fun in (:vev_KI, :vev_IK, :tr)
        Ovac = I_H
        p0 = model_val(Ovac, Ovac, fun)
        # P1 全行核对（K=真空：O_K = I）
        p1 = [model_val(I_H, cascade_img([(A, 1.0 + 0.0im)]; order=ord), fun) for A in EVENS]
        e1 = maximum(abs.(p1 - row1))
        p2 = [model_val(cascade_img([(A, 1.0 + 0.0im)]; order=ord), I_H, fun) for A in EVENS]
        e2 = maximum(abs.(p2 - row2))
        @printf("bar→c†=%-5s order=:%-5s fun=:%-6s : P0=%.3f  P1err=%.2e  P2err=%.2e\n",
                bar2c, string(ord), string(fun), p0, e1, e2)
    end

    # ------------------------------------------------ P3：εd 依赖的单门 --
    println("\n== P3：单门级联 vs 锚（K 链门 c=e^{−τε}−1）==")
    ring = [0.5, 0.5]     # order-2 均匀 ring
    for ed in [0.0, 0.5, 1.0, 2.0]
        cg = exp(-ring[1] * ed) - 1
        # K 链：与 mk_kchain 同构（order-2 nl=2：gates @slice 2,3 + resolutions + BC）
        factors = Tuple{Vector{Int},ComplexF64}[]
        for p in 1:nl
            iszero(ed) || push!(factors, (sort([spos[(p+1, true)], spos[(p+1, false)]]),
                                          ComplexF64(exp(-ring[p] * ed) - 1)))
            s1 = p + 1; s2 = p == nl ? 2 : p + 2
            push!(factors, (sort([spos[(s1, true)], spos[(s2, false)]]), -1.0 + 0.0im))
        end
        push!(factors, (sort([spos[(1, true)], spos[(k, false)]]), -1.0 + 0.0im))
        # GTEMPO 侧（同 mk_kchain 序）
        mpsK = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (A, c) in factors
            mpsK = GTEMPO.apply!(exp(GTEMPO.GTerm([p + 1 for p in A]...; coeff=ComplexF64(c))), mpsK)
        end
        m = [meas(mpsK, inj_mps[A]) for A in EVENS[1:6]]
        println("εd=$ed 实测 [B=∅,单对×5]：", m)
    end
    println("\n（模型侧对 P3 的复算需先由 P0-P2 定出约定——见上表）")
end
main()
