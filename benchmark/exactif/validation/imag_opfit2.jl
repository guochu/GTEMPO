# =============================================================================
# imag_opfit2.jl —— 位置模代数模型判决（nl=2 穷举，参数正确的版本）
#
# 模型：生成元 θ_p ↔ 位置 p 的费米模式 a_p（CAR，2L 模式，dim 2^{2L}=64）。
#   因子 exp(GTerm(A,c)) ↔ I + c·∏_{p∈A} a_p（升序积）；
#   级联 ↔ 代数有序积 M；复合 integrate(K,I) ↔ Λ(M_K·M_I)（或反序），
#   Λ = 固定线性泛函，未知量 = Λ 在 32 个偶单项式上的值。
# 判决：
#   1. τ2 表（1024 方程）解 32 个 Λ_b —— 真实一致性检验
#   2. 解出的 Λ 预测级联探针（K(0) 链 + 随机级联）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_opfit2.jl
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

function position_modes(L::Int)
    dim = 1 << L
    am = [zeros(ComplexF64, dim, dim) for _ in 1:L]
    for n in 0:(dim-1), p in 1:L
        ((n >> (p-1)) & 1) == 1 || continue
        jw = isodd(count_ones(n & ((1 << (p-1)) - 1))) ? -1.0 : 1.0
        am[p][(n ⊻ (1 << (p-1))) + 1, n + 1] = jw
    end
    return am
end


"""朴素循环矩阵乘（绕开 OpenBLAS zgemm 在本机的自旋 bug）。"""
function muln(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    da, db = size(A, 1), size(B, 2)
    dk = size(A, 2)
    C = zeros(ComplexF64, da, db)
    @inbounds for i in 1:da, k in 1:dk
        aik = A[i, k]
        aik == 0 && continue
        @simd for j in 1:db
            C[i, j] += aik * B[k, j]
        end
    end
    return C
end

function main()
    nl = 2
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT); np = 2L; dim = 1 << np
    am = position_modes(np)
    println(stderr, "[t] modes done"); flush(stderr)
    I_M = Matrix{ComplexF64}(I, dim, dim)

    factor_img = (A::Vector{Int}, c::ComplexF64) -> begin
        X = I_M * 0.0; first = true
        for p in sort(A)
            X = first ? am[p+1] : muln(X, am[p+1])
            first = false
        end
        return I_M + c .* X
    end
    cascade_img = function(fs::Vector{Tuple{Vector{Int},ComplexF64}})
        O = I_M
        for (A, c) in fs
            O = muln(factor_img(A, c), O)
        end
        return O
    end

    EVENS = Vector{Vector{Int}}([Int[]])
    for m in 2:2:L, comb in combinations(0:(L-1), m)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    println(stderr, "[t] evens ", NE); flush(stderr)

    # 偶单项式算符 a_b 与 CAR 对偶 D_b
    mono_op = function(A::Vector{Int})
        isempty(A) && return I_M
        X = I_M * 0.0; first = true
        for p in sort(A)
            X = first ? am[p+1] : muln(X, am[p+1])
            first = false
        end
        return X
    end
    @time Aops = [mono_op(A) for A in EVENS]
    Dops = [Matrix(adjoint(Aops[i])) for i in 1:NE]
    @time G = [tr(muln(Dops[i], Aops[j])) for i in 1:NE, j in 1:NE]
    println("对偶度规 G 对角性：max offdiag = ", round(maximum(abs.(G .- Diagonal(diag(G)))), digits=6),
            "   diag 范围 = ", round.(abs.(diag(G)), digits=3))

    # m_b(M) = Tr(D_b M)/G[b,b] —— 多项式系数提取
    mcoef = function(M::Matrix{ComplexF64})
        return [tr(muln(Dops[i], M)) / G[i, i] for i in 1:NE]
    end

    # ------------------------------------------------ 载入 τ2 表 --
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    tau2 = Array{ComplexF64,2}(undef, NE, NE)
    for (i, ln) in enumerate(eachline(joinpath(outdir, "imag_tau2.csv")))
        i == 1 && continue
        f = split(strip(ln), ",")
        tau2[parse(Int, f[1]), parse(Int, f[2])] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
    end
    println("τ2 表载入完成（|τ2| 非零率 = ", round(count(x -> abs(x) > 1e-9, tau2) / NE^2, digits=3), "）")

    # 级联像：M_a = ∏(I + a_pair)（injected 序）
    @time Ma = [cascade_img([(collect(pair), 1.0 + 0.0im) for pair in
                       [(A[2i-1], A[2i]) for i in 1:(length(A)÷2)]]) for A in EVENS]
    println(stderr, "[t] Ma done"); flush(stderr)

    # -------------------------------------- 步骤 1：32 未知 Λ_b 拟合 --
    Fm = Matrix{ComplexF64}(undef, NE * NE, NE)
    z = Vector{ComplexF64}(undef, NE * NE)
    li = 0
    for ia in 1:NE, ib in 1:NE
        li += 1
        Fm[li, :] .= mcoef(muln(Ma[ia], Ma[ib]))
        z[li] = tau2[ia, ib]
    end
    Λb = (Fm' * Fm + 1e-12 * I) \ (Fm' * z)
    res = norm(Fm * Λb - z) / max(norm(z), 1e-30)
    @printf("Λ 拟合相对残差 = %.3e（32 未知 vs %d 方程）%s\n", res, NE * NE,
            res < 1e-8 ? "✅ 一致" : "✗ 不一致")
    println("Λ_b（偶单项式泛函值，非零项）：")
    for (i, A) in enumerate(EVENS)
        abs(Λb[i]) > 1e-9 || continue
        @printf("   Λ(%s) = %8.4f%+8.4fim\n", join(A, ","), real(Λb[i]), imag(Λb[i]))
    end

    # --------------------------------- 步骤 2：级联探针 holdout --
    spos = Dict{Tuple{Int,Bool}, Int}()
    for s in 1:nl+1
        spos[(s, true)] = LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1
        spos[(s, false)] = LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1
    end
    k0 = Tuple{Vector{Int},ComplexF64}[]
    for p in 1:nl
        s1 = p + 1; s2 = p == nl ? 2 : p + 2
        push!(k0, (sort([spos[(s1, true)], spos[(s2, false)]]), -1.0 + 0.0im))
    end
    push!(k0, (sort([spos[(1, true)], spos[(nl + 1, false)]]), -1.0 + 0.0im))

    inj = function(A::Vector{Int})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for i in 1:2:length(A)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(A[i] + 1, A[i+1] + 1; coeff=ComplexF64(1.0))), mps)
        end
        return mps
    end
    vac = GTEMPO.vacuumstate(ComplexF64, LAT)
    meas = function(fs)
        mps = vac
        for (A, c) in fs
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm([p + 1 for p in A]...; coeff=ComplexF64(c))), mps)
        end
        return mps
    end
    Λpred = function(fsK, fsI; order::Symbol=:KI)
        MK = cascade_img(fsK); MI = cascade_img(fsI)
        M = order === :KI ? muln(MK, MI) : muln(MI, MK)
        return dot(mcoef(M), Λb)
    end

    println("\n级联探针判决：")
    z0 = real(GTEMPO.integrate(LAT, meas(k0), vac))
    @printf("K(0)+vac：实测 %.6f   模型 KI %.6f   模型 IK %.6f\n",
            z0, real(Λpred(k0, Tuple{Vector{Int},ComplexF64}[]; order=:KI)),
            real(Λpred(k0, Tuple{Vector{Int},ComplexF64}[]; order=:IK)))

    rng = MersenneTwister(7)
    allpairs = [(a, b) for b in 1:L for a in 1:(b-1)]
    errKI = 0.0; errIK = 0.0
    for _ in 1:30
        fsK = Tuple{Vector{Int},ComplexF64}[]
        for _ in 1:(rand(rng, 1:3))
            (a, b) = rand(rng, allpairs)
            push!(fsK, (sort([a - 1, b - 1]), ComplexF64(randn(rng))))
        end
        B = EVENS[rand(rng, 1:NE)]
        fsI = [(collect(pair), 1.0 + 0.0im) for pair in [(B[2i-1], B[2i]) for i in 1:(length(B)÷2)]]
        zr = GTEMPO.integrate(LAT, meas(fsK), inj(B))
        errKI = max(errKI, abs(real(Λpred(fsK, fsI; order=:KI)) - real(zr)))
        errIK = max(errIK, abs(real(Λpred(fsK, fsI; order=:IK)) - real(zr)))
    end
    @printf("随机级联 30 探针：maxerr KI = %.2e   IK = %.2e\n", errKI, errIK)
    ok = (errKI < 1e-9 || errIK < 1e-9) && res < 1e-8
    println(ok ? "✅ 位置模代数模型成立（可升 nl=6）" : "✗ 模型不成立（需更深结构）")
end
main()
