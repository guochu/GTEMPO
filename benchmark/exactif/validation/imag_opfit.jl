# =============================================================================
# imag_opfit.jl —— φ 双线性配对反推（Plan C，nl=2 小格）
#
# 假设：integrate(K, I) = φ(O_K, O_I) = Σ_{αβγδ} T[αβγδ]·OK[αβ]·OI[γδ]
#   O = 因子级联的有序积像；嵌入方向 × 级联序为待扫描约定。
# 方法：随机级联探针 × 全 32 个 I 基矢 → 4 指标设计矩阵线性最小二乘解 T
#   （dim⁴ = 4096 复参数）→ holdout 验证 → 与解析泛函候选对照。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_opfit.jl
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0

function mode_ops(k::Int)
    dim = 1 << k
    cs = [zeros(ComplexF64, dim, dim) for _ in 1:k]
    for n in 0:(dim-1), s in 1:k
        ((n >> (s-1)) & 1) == 1 || continue
        jw = isodd(count_ones(n & ((1 << (s-1)) - 1))) ? -1.0 : 1.0
        cs[s][(n ⊻ (1 << (s-1))) + 1, n + 1] = jw
    end
    return cs, [adjoint(cs[s]) for s in 1:k]
end

"""4 指标特征：f[α,β,γ,δ] = A[α,β]·B[γ,δ] → vec。"""
function feat4(A::Matrix{ComplexF64}, B::Matrix{ComplexF64}, dim::Int)
    A4 = reshape(A, dim, 1, dim, 1)
    B4 = reshape(B, 1, dim, 1, dim)
    return vec(permutedims(A4 .* B4, (1, 3, 2, 4)))
end

function main()
    nl = 2
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT); k = nl + 1; dim = 1 << k
    spos = Dict{Tuple{Int,Bool}, Int}()
    for s in 1:k
        spos[(s, true)] = LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1
        spos[(s, false)] = LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1
    end
    pos2leg = Dict(v => kk for (kk, v) in spos)
    cs, cds = mode_ops(k)
    I_H = Matrix{ComplexF64}(I, dim, dim)
    vac0 = zeros(ComplexF64, dim); vac0[1] = 1.0

    factor_img = function(A::Vector{Int}, c::ComplexF64; bar2c::Bool)
        X = I_H * 0.0; first = true
        for p in sort(A)
            (sl, bar) = pos2leg[p]
            Xp = bar ? (bar2c ? cds[sl] : cs[sl]) : (bar2c ? cs[sl] : cds[sl])
            X = first ? Xp : X * Xp
            first = false
        end
        return I_H + c .* X
    end
    cascade_img = function(fs; bar2c, order)
        O = I_H
        for (A, c) in fs
            F = factor_img(A, c; bar2c=bar2c)
            O = order === :left ? F * O : O * F
        end
        return O
    end

    # ------------------------------------------------------------ GTEMPO 探针 --
    EVENS = Vector{Vector{Int}}([Int[]])
    for m in 2:2:L, comb in combinations(0:(L-1), m)
        push!(EVENS, collect(comb))
    end
    inj = function(A::Vector{Int})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for i in 1:2:length(A)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(A[i] + 1, A[i+1] + 1; coeff=ComplexF64(1.0))), mps)
        end
        return mps
    end
    inj_mps = [inj(A) for A in EVENS]
    vac = GTEMPO.vacuumstate(ComplexF64, LAT)

    rng = MersenneTwister(20261108)
    allpairs = [(a, b) for b in 1:L for a in 1:(b-1)]
    Klists = Vector{Vector{Tuple{Vector{Int},ComplexF64}}}()
    for _ in 1:100
        nf = rand(rng, 1:3)
        fs = Vector{Tuple{Vector{Int},ComplexF64}}()
        for _ in 1:nf
            (a, b) = rand(rng, allpairs)
            push!(fs, (sort([a - 1, b - 1]), ComplexF64(randn(rng) + im * randn(rng))))
        end
        push!(Klists, fs)
    end
    # K(0) 原型（resolution+BC，含 O(1) 系数多样性变体）
    k0 = Tuple{Vector{Int},ComplexF64}[]
    for p in 1:nl
        s1 = p + 1; s2 = p == nl ? 2 : p + 2
        push!(k0, (sort([spos[(s1, true)], spos[(s2, false)]]), -1.0 + 0.0im))
    end
    push!(k0, (sort([spos[(1, true)], spos[(k, false)]]), -1.0 + 0.0im))
    push!(Klists, k0)

    meas = function(fs)
        mps = vac
        for (A, c) in fs
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm([p + 1 for p in A]...; coeff=ComplexF64(c))), mps)
        end
        return mps
    end

    println("测量探针复合量（$(length(Klists)) 列表 × $(length(EVENS)) 基矢）...")
    Zval = Vector{Vector{ComplexF64}}()
    mpsKs = Vector{Any}()
    for fs in Klists
        mpsK = meas(fs)
        push!(mpsKs, mpsK)
        push!(Zval, [GTEMPO.integrate(LAT, mpsK, injm) for injm in inj_mps])
    end
    println("测量完成")

    perm = randperm(rng, length(Klists))
    ntrain = length(Klists) - 15
    train_idx = perm[1:ntrain]; hold_idx = perm[ntrain+1:end]

    for bar2c in (true, false), ord in (:left, :right)
        OI_mats = [cascade_img([(collect(pair), 1.0 + 0.0im) for pair in
                                [(A[2i-1], A[2i]) for i in 1:(length(A)÷2)]]; bar2c=bar2c, order=ord)
                   for A in EVENS]
        OK_mats = [cascade_img(Klists[i]; bar2c=bar2c, order=ord) for i in eachindex(Klists)]

        F_train = Vector{Vector{ComplexF64}}(); z_train = ComplexF64[]
        for ki in train_idx, (bi, A) in enumerate(EVENS)
            push!(F_train, feat4(OK_mats[ki], OI_mats[bi], dim))
            push!(z_train, Zval[ki][bi])
        end
        M = permutedims(hcat(F_train...))
        Tfit = (M' * M + 1e-10 * I) \ (M' * z_train)
        res = norm(M * Tfit - z_train) / max(norm(z_train), 1e-30)
        eh = 0.0
        for ki in hold_idx, (bi, A) in enumerate(EVENS)
            pred = dot(feat4(OK_mats[ki], OI_mats[bi], dim), Tfit)
            eh = max(eh, abs(pred - Zval[ki][bi]))
        end
        @printf("bar2c=%-5s ord=:%-5s : train relres = %.2e   holdout maxerr = %.2e\n",
                bar2c, string(ord), res, eh)

        if res < 1e-8
            Tm = reshape(Tfit, dim, dim, dim, dim)
            nz = count(x -> abs(x) > 1e-8, Tm)
            @printf("   T 非零元（>1e-8）= %d / %d\n", nz, length(Tm))
            # 解析候选对照
            cands = Dict{String,Array{ComplexF64,4}}()
            c1 = zeros(ComplexF64, dim, dim, dim, dim)     # Tr(A·B)
            for μ in 1:dim, ν in 1:dim; c1[μ, ν, ν, μ] = 1.0; end
            cands["Tr(A·B)"] = c1
            c2 = zeros(ComplexF64, dim, dim, dim, dim)     # ⟨0|A·B|0⟩
            for ν in 1:dim; c2[1, ν, ν, 1] = 1.0; end
            cands["⟨0|A·B|0⟩"] = c2
            c3 = zeros(ComplexF64, dim, dim, dim, dim)     # ⟨0|B·A|0⟩
            for ν in 1:dim; c3[ν, 1, 1, ν] = 1.0; end
            cands["⟨0|B·A|0⟩"] = c3
            for (nm, Tc) in cands
                @printf("   ‖T − %s‖∞ = %.2e\n", nm, maximum(abs.(Tm .- Tc)))
            end
            # T 的非零模式（行/列指标分布）
            idxs = findall(x -> abs(x) > 1e-8, Tm)
            println("   非零 T 元样本（前 20）：")
            for t in idxs[1:min(20, end)]
                @printf("     T[%d,%d,%d,%d] = %8.3f%+8.3fim\n", Tuple(t)..., real(Tm[t]), imag(Tm[t]))
            end
        end
    end
    # K(0)+vac 判据
    println("\nK(0)+vac 实测 = ", real(GTEMPO.integrate(LAT, meas(k0), vac)),
            "（各约定下 ⟨0|O_K(0)|0⟩ / Tr：见上——若全不合则 φ 需非平凡结构）")
end
main()
