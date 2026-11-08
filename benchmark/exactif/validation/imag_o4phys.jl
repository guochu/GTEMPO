# =============================================================================
# imag_o4phys.jl —— 物理子空间限制拟合（盲拟合的最后一搏/判决）
#
# 原理：行空间的基方向 u_j 的采样函数 f_j = A'u_j 可按"核单项式族"分类——
#   物理：f_j ∈ span{x_a^p x_c^q, p≤2, q≤4}（门乘积代数，15 维，无增长方向）
#   伪：f_j 需要有理/增长分量（残差大）
# 限制在物理方向上解约束 → Ist 系数量级恢复 → 网格外行为判决。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4phys.jl
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
    F = svd(A)

    # 核单项式族设计矩阵：x_a = e^{−0.6756εd}, x_c = e^{+0.0878εd}（τ 取负注意符号）
    # ring 时长：p 门 τ_a>0（×2），c 门 τ_c<0（×4）
    ta = ring_links[1]; tc = ring_links[2]
    @printf("ring 时长：τ_a=%.4f (×%d), τ_c=%.4f (×%d)\n", ta, count(==(ta), ring_links), tc, count(==(tc), ring_links))
    monos = [(p, q) for p in 0:2 for q in 0:4]
    Φ = [ComplexF64(exp(-ta * ed)^p * exp(-tc * ed)^q) for ed in vcat(heds, geds), (p, q) in monos]

    println("行空间基方向分类（f_j = A'u_j 对 15 维核单项式族的拟合残差）：")
    physP = Int[]; quirkP = Int[]
    for j in eachindex(F.S)
        σj = F.S[j]
        σj > 1e-14 || continue
        fj = A' * F.U[:, j]                # 36 向量（采样函数）
        c = Φ \ fj
        rj = norm(fj - Φ * c) / max(norm(fj), 1e-30)
        if rj < 1e-9
            push!(physP, j)
        else
            push!(quirkP, j)
            @printf("   伪方向 σ%-3d = %10.3e : 相对残差 %8.2e   |f| = %8.2e\n", j, σj, rj, norm(fj))
        end
    end
    @printf("物理方向 %d 个 / 伪方向 %d 个\n", length(physP), length(quirkP)); flush(stdout)

    # 物理子空间限制解
    Ist = zeros(ComplexF64, NE)
    for j in physP
        Ist .+= (dot(F.V[:, j], rhs) / F.S[j]) .* F.U[:, j]
    end
    r = A' * Ist - rhs
    @printf("物理子空间解：残差 = %.3e   ‖Ist‖∞ = %.3e\n", norm(r), maximum(abs, Ist))
    println("最大逐行残差：")
    perm = sortperm(abs.(r), rev=true)
    for i in perm[1:min(6, count(!iszero, r))]
        abs(r[i]) > 1e-12 && @printf("   %-12s : %10.3e\n", tags[i], abs(r[i]))
    end

    # rhs 在伪方向上的分量（若大 → Z 需要伪方向，物理子空间必失败）
    vr = F.V' * rhs
    rq = sum(abs2, [vr[j] for j in quirkP])
    @printf("rhs 在伪方向分量平方和的平方根 = %.3e\n", sqrt(rq))

    # 离线导出 Ist 供 MPS 构建（物理解，ν 归一在 closure 内做）
    open(joinpath(outdir, "imag_n2o4ma_ist_phys.txt"), "w") do io
        for (i, Aset) in enumerate(EVENS)
            abs(Ist[i]) > 1e-12 || continue
            println(io, join(sort(Aset), " "), " ", real(Ist[i]), " ", imag(Ist[i]))
        end
    end
    println("物理解已导出 imag_n2o4ma_ist_phys.txt")
end
main()
