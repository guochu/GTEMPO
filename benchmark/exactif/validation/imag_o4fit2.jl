# =============================================================================
# imag_o4fit2.jl —— order-4 "不一致"是结构缺口还是正则化伪影？
#
# 假设：非均匀调度的 ring 门含小时长（|c·δt| ~ 0.09），(e^{−τε}−1) 内核系数
# ~1e-2，4 个 c 门乘积 ~1e-5 → h 行矩阵病态；closure 的 ridge=1e-13 加在
# 病态 A'A 上淹没小方向 → 假不一致。
#
# 判决：SVD 谱 + 截断 pinv 拟合残差 vs ridge 拟合残差。
#   截断 pinv 把残差打到 ~1e-12 → 伪影（改拟合器即闭环）
#   截断 pinv 仍 ~2e-5      → 真结构缺口
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4fit2.jl [tag...]
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
    for part in (1, 2)
        for (i, ln) in enumerate(eachline(joinpath(outdir, "imag_$(tag)_h$part.csv")))
            i == 1 && continue
            f = split(strip(ln), ",")
            get!(hrows, parse(Float64, f[1]), zeros(ComplexF64, NE))[parse(Int, f[2])] =
                ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
        end
    end
    grows = Dict{Float64, Vector{ComplexF64}}()
    for part in (1, 2)
        for (i, ln) in enumerate(eachline(joinpath(outdir, "imag_$(tag)_g$part.csv")))
            i == 1 && continue
            f = split(strip(ln), ",")
            g = get!(grows, parse(Float64, f[1]), zeros(ComplexF64, 3NE))
            ib = parse(Int, f[2])
            g[ib] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
            g[NE + ib] = ComplexF64(parse(Float64, f[5]), parse(Float64, f[6]))
            g[2NE + ib] = ComplexF64(parse(Float64, f[7]), parse(Float64, f[8]))
        end
    end
    return hrows, grows
end

function run_tag(outdir, tag, Nmacro, order)
    # env 覆盖调度（isobatch 数据诊断用：IMAG_CELLS/IMAG_RLINKS）
    if haskey(ENV, "IMAG_CELLS") && haskey(ENV, "IMAG_RLINKS")
        cells_env = [parse(Float64, x) for x in split(ENV["IMAG_CELLS"], ",")]
        rlinks_env = [parse(Float64, x) for x in split(ENV["IMAG_RLINKS"], ",")]
        sched = (cells = cells_env, links = rlinks_env)
        ring_links = rlinks_env
    else
        sched = NativeExact.native_schedule(order, Nmacro, BETA / Nmacro)
        ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    end
    nl = length(ring_links)
    L = 2 * (nl + 1)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    println("="^72)
    @printf("tag=%s  N=%d order=%d  nl=%d L=%d evensets=%d\n", tag, Nmacro, order, nl, L, NE)
    flush(stdout)

    hrows, grows = load_rows(outdir, tag, NE)
    heds = sort(collect(keys(hrows)))
    geds = sort(collect(keys(grows)))
    cells = sched.cells

    rowsH = [hrows[ed] for ed in heds]
    rhsH = [fock_anchor(cells, ring_links, ed)[1] for ed in heds]
    rowsG = [hrows[ed] .- grows[ed][NE .+ (1:NE)] .- real(fock_anchor(cells, ring_links, ed)[2]) .* hrows[ed] for ed in geds]
    A = hcat(vcat(rowsH, rowsG)...)
    rhs = vcat(rhsH, fill(0.0 + 0.0im, length(geds)))
    @printf("A: %d×%d\n", size(A, 1), size(A, 2)); flush(stdout)

    # 1) SVD 谱
    F = svd(A)
    sv = F.S
    @printf("σ 谱：σ1=%.3e σ2=%.3e σ3=%.3e  …  σ12=%.3e σ13=%.3e σ14=%.3e … σmin=%.3e\n",
            sv[1], sv[2], sv[3], sv[min(12, end)], sv[min(13, end)], sv[min(14, end)], sv[end])
    println("  全谱（科学计数）：")
    for i in eachindex(sv)
        @printf("   σ%-3d = %10.3e%s", i, sv[i], i % 4 == 0 ? "\n" : "")
    end
    println()

    # 2) ridge 拟合（closure 现行）
    Ist_r = A * ((A' * A + 1e-13 * I) \ rhs)
    println(@sprintf("ridge 1e-13 残差 = %.3e", norm(A' * Ist_r - rhs)))

    # 3) 截断 SVD 拟合：Ist = Σ_{σ>tol} (v'rhs)/σ · u
    println("\n截断 SVD 拟合（min‖A'I−rhs‖）：")
    for tol in (1e-16, 1e-14, 1e-12, 1e-10, 1e-8)
        keep = F.S .> tol
        Ist = F.U * ((F.V' * rhs) ./ ifelse.(keep, F.S, Inf))
        res = norm(A' * Ist - rhs)
        nkeep = count(keep)
        @printf("tol=%.0e : keep=%2d  残差 = %.3e\n", tol, nkeep, res)
    end
    # 真距离：rhs 到 col(A')（=rhs 可被模型预测的部分）
    proj = F.V * (F.V' * rhs)             # rhs 在 row-space（v 向量）上的投影
    @printf("rhs 到 row-space(A) 真距离 = %.3e   (|rhs| = %.3e)\n", norm(rhs - proj), norm(rhs))
    # 残差分解到小 σ 方向
    vr = F.V' * rhs
    println("rhs 在各右奇异方向上的分量 |⟨v_i, rhs⟩|：")
    for i in eachindex(F.S)
        @printf("   σ%-3d = %10.3e : |comp| = %10.3e%s", i, F.S[i], abs(vr[i]), i % 3 == 0 ? "\n" : "")
    end
    println()
    flush(stdout)
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    tags = isempty(ARGS) ? ["n1o4ma", "n2o4ma"] : ARGS
    specs = Dict("n2o2ma" => (2, 2), "n2o1ma" => (2, 1), "n3o2ma" => (3, 2),
                 "n1o4ma" => (1, 4), "n2o4ma" => (2, 4), "n1o2ma" => (1, 2))
    for tag in tags
        (N, o) = haskey(ENV, "IMAG_CELLS") ? (1, 2) : specs[tag]
        run_tag(outdir, tag, N, o)
    end
end
main()
