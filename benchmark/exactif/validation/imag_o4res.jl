# =============================================================================
# imag_o4res.jl —— 组合锚点（旧 24 + 大 εd 12）残差解剖
#
# 问题：36 约束系统 keep=34 残差 1.8e-6（两个 σ≤1e-14 方向）。定位：
#   1. σ 谱全表（找出被截断的方向索引）
#   2. 逐锚点残差表（残差集中在旧点还是新点、h 还是 g）
#   3. 被截断方向的 u 向量在锚点上的分布
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4res.jl
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
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

annih_op_f(dim::Int, mode::Int) = begin
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    op
end
function fock_anchor(cells, gate_durs, eps_d)
    m = length(VB); nt = 1 + m; dim = 1 << nt
    ops = [annih_op_f(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim)
    h = hcoupled()
    for i in 1:nt, j in 1:nt
        iszero(h[i, j]) || (Hq .+= h[i, j] .* (ops[i]' * ops[j]))
    end
    Him = eps_d .* (ops[1]' * ops[1])
    Z = Matrix{ComplexF64}(I, dim, dim); Zn = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z; Zn = W * Zn
    end
    Zn = ops[1]' * ops[1] * Zn
    return tr(Z), tr(Zn) / tr(Z)
end

function load_rows(outdir, tag, NE)
    hrows = Dict{Float64, Vector{ComplexF64}}()
    for file in sort(filter(f -> occursin("imag_$(tag)_h", f) && endswith(f, ".csv"), readdir(outdir, join=true)))
        for (i, ln) in enumerate(eachline(file))
            i == 1 && continue
            f = split(strip(ln), ",")
            get!(hrows, parse(Float64, f[1]), zeros(ComplexF64, NE))[parse(Int, f[2])] =
                ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
        end
    end
    grows = Dict{Float64, Vector{ComplexF64}}()
    for file in sort(filter(f -> occursin("imag_$(tag)_g", f) && endswith(f, ".csv"), readdir(outdir, join=true)))
        for (i, ln) in enumerate(eachline(file))
            i == 1 && continue
            f = split(strip(ln), ",")
            ed = parse(Float64, f[1]); ib = parse(Int, f[2])
            g = get!(grows, ed, zeros(ComplexF64, 3NE))
            g[ib] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
            g[NE + ib] = ComplexF64(parse(Float64, f[5]), parse(Float64, f[6]))
            g[2NE + ib] = ComplexF64(parse(Float64, f[7]), parse(Float64, f[8]))
        end
    end
    return hrows, grows
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    tag = "n2o4ma"
    NL, order = 2, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    L = 2 * (nl + 1)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    hrows, grows = load_rows(outdir, tag, NE)
    heds = sort(collect(keys(hrows)))
    geds = sort(collect(keys(grows)))
    cells = sched.cells

    tags = String[]; rowsV = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]
    for ed in heds
        push!(rowsV, hrows[ed]); push!(rhs, fock_anchor(cells, ring_links, ed)[1]); push!(tags, "Z@$ed")
    end
    for ed in geds
        push!(rowsV, hrows[ed] .- grows[ed][NE .+ (1:NE)] .- real(fock_anchor(cells, ring_links, ed)[2]) .* hrows[ed])
        push!(rhs, 0.0 + 0.0im); push!(tags, "n@$ed")
    end
    A = hcat(rowsV...)
    @printf("A: %d×%d\n", size(A, 1), size(A, 2))
    F = svd(A)
    println("σ 谱（全部 36）：")
    for i in eachindex(F.S)
        @printf("   σ%-3d = %10.3e%s", i, F.S[i], i % 4 == 0 ? "\n" : "")
    end
    println()

    Ist = F.U * ((F.V' * rhs) ./ ifelse.(F.S .> 1e-14, F.S, Inf))
    r = A' * Ist - rhs
    @printf("keep>1e-14 残差 = %.3e\n", norm(r))
    println("逐锚点残差：")
    for i in eachindex(tags)
        abs(r[i]) > 1e-9 && @printf("   %-12s : %10.3e\n", tags[i], abs(r[i]))
    end

    # rhs 在被截断方向上的分量
    vr = F.V' * rhs
    println("rhs 分量 vs σ（|⟨v_i,rhs⟩| > 1e-9 或 σ < 1e-13）：")
    for i in eachindex(F.S)
        (F.S[i] < 1e-13 || abs(vr[i]) > 1e-9) &&
            @printf("   σ%-3d = %10.3e : |comp| = %10.3e\n", i, F.S[i], abs(vr[i]))
    end

    # 被截断方向 u 的锚点分布（|u_i| 大的锚点）
    for i in findall(F.S .<= 1e-14)
        println("截断方向 σ$i 的 |u| 分布（>0.15 的锚点）：")
        for j in eachindex(tags)
            abs(F.V[j, i]) > 0.15 && @printf("      %-12s : %8.3f\n", tags[j], abs(F.V[j, i]))
        end
    end

    # 新旧行的自洽交叉检验：大 εd 行是否与旧行的函数空间一致
    # 用旧 24 行的 SVD 拟合新 12 行的 rhs
    olds = [i for (i, tg) in enumerate(tags) if !occursin("3.0", tg) && !occursin("4.5", tg) &&
                                             !occursin("6.0", tg) && !occursin("8.0", tg) &&
                                             !occursin("10.0", tg) && !occursin("12.0", tg)]
    Aold = A[:, olds]
    Fold = svd(Aold)
    rhsold = rhs[olds]; rhsnew = rhs[setdiff(eachindex(tags), olds)]
    Istold = Fold.U * ((Fold.V' * rhsold) ./ ifelse.(Fold.S .> 1e-14, Fold.S, Inf))
    projnew = A[:, setdiff(eachindex(tags), olds)]' * Istold
    println("\n旧 24 行拟合解对新 12 行的预测残差：")
    news = setdiff(eachindex(tags), olds)
    for (j, i) in enumerate(news)
        abs(projnew[j] - rhsnew[j]) > 1e-9 &&
            @printf("   %-12s : 预测 %10.6f vs 实际 %10.6f (差 %10.3e)\n",
                    tags[i], real(projnew[j]), real(rhsnew[j]), abs(projnew[j] - rhsnew[j]))
    end
    @printf("新行预测残差范数 = %.3e\n", norm(projnew - rhsnew))
end
main()
