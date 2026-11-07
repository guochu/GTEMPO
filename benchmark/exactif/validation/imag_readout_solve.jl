# =============================================================================
# imag_readout_solve.jl —— n-读出关系发现：联合一致性最小化
#
# 事实：全 32 系数一般 I 也无法同时满足 Z-行 + {1−g_end} n-行（残差 1.8e-6）
# ⇒ 现行"n = 1 − g_end"读出公式错误。本脚本：n_pipe = α·g_end + β·g_1（γ=0）
# 的 (α,β) 复二维一致性最小化（含网格 + 精化），并报告一致解下的 I 结构。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_readout_solve.jl
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))

const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = begin
    m = length(VB)
    h = zeros(ComplexF64, m + 1, m + 1)
    for k in 1:m
        h[1, k+1] = conj(VB[k]); h[k+1, 1] = VB[k]; h[k+1, k+1] = HB[k, k]
    end
    h
end
annih_op(dim::Int, mode::Int) = begin
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
    ops = [annih_op(dim, k) for k in 0:(nt-1)]
    Hq = zeros(ComplexF64, dim, dim); h = hcoupled()
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

const NL = 2
const N = 2
const CELLS = fill(BETA / N, N)
const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)

function load_csv(path; cols::Int)
    out = Dict{Tuple{Int,Float64}, Matrix{ComplexF64}}()
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        f = split(strip(ln), ",")
        key = (parse(Int, f[1]), parse(Float64, f[2]))
        ib = parse(Int, f[3])
        M = get!(out, key, zeros(ComplexF64, 32, cols))
        for c in 1:cols
            M[ib, c] = ComplexF64(parse(Float64, f[2 + 2c]), parse(Float64, f[3 + 2c]))
        end
    end
    return out
end

function main()
    hdata = load_csv(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_h.csv"); cols=1)
    gdata = load_csv(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_g.csv"); cols=3)
    eds = sort([k[2] for k in keys(gdata) if k[1] == 2])
    hed_all = sort([k[2] for k in keys(hdata) if k[1] == 2])

    """联合系统残差（给定读出系数 α, β）：[Z-行; (α g_end + β g_1 − n·z)-行]。"""
    function joint_min_residual(alpha::ComplexF64, beta::ComplexF64)
        rows = Vector{ComplexF64}[]; rhs = ComplexF64[]
        for ed in hed_all
            push!(rows, hdata[(2, ed)][:, 1]); push!(rhs, fock_anchor(CELLS, CELLS, ed)[1])
        end
        for (r, ed) in enumerate(eds)
            hz = hdata[(2, ed)][:, 1]
            gv = gdata[(2, ed)]
            nanc = real(fock_anchor(CELLS, CELLS, ed)[2])
            nrow = alpha .* gv[:, 2] .+ beta .* gv[:, 3] .- nanc .* hz
            push!(rows, nrow); push!(rhs, 0.0 + 0.0im)
        end
        A = hcat(rows...)
        sol = A * ((A' * A + 1e-12 * I) \ b_rhs(rhs))
        return norm(A' * sol - b_rhs(rhs)), sol
    end
    b_rhs(r) = r

    # (α, β) 网格（实部/虚部 ∈ [-2, 2]）粗搜 → 最优点精化
    best = (Inf, 0.0 + 0.0im, 0.0 + 0.0im, nothing)
    for ar in -2.0:0.25:2.0, ai in -1.0:0.25:1.0, br in -2.0:0.25:2.0, bi in -1.0:0.25:1.0
        rf, _ = joint_min_residual(ComplexF64(ar, ai), ComplexF64(br, bi))
        rf < best[1] && (best = (rf, ComplexF64(ar, ai), ComplexF64(br, bi), nothing))
    end
    println("粗网格最优 |F| = ", best[1], "  α = ", best[2], "  β = ", best[3])
    # 局部精化（随机局部搜索）
    Random.seed!(3)
    for it in 1:400
        a2 = best[2] + ComplexF64(0.02randn(), 0.02randn())
        b2 = best[3] + ComplexF64(0.02randn(), 0.02randn())
        rf, _ = joint_min_residual(a2, b2)
        if rf < best[1]
            best = (rf, a2, b2, nothing)
        end
    end
    println("精化最优 |F| = ", best[1], "  α = ", best[2], "  β = ", best[3])
    if best[1] < 1e-8
        rf, sol = joint_min_residual(best[2], best[3])
        println("\n一致解找到：|F| = $rf")
        println("n-读出公式：n_pipe = α·g_end + β·g_1,  α = ", best[2], "  β = ", best[3])
        println("\n对应 min-norm I 的非零项：")
        for (i, Aset) in enumerate(EVENS_RE)
            v = sol[i]
            abs(v) > 1e-8 || continue
            @printf("   %-16s : %11.7f%s\n", Aset, real(v),
                    abs(imag(v)) > 1e-8 ? " + $(round(imag(v), digits=7))im" : "")
        end
    end
end
const EVENS_RE = begin
    L = 6
    ev = Vector{Vector{Int}}([Int[]])
    for k in (2, 4, 6), comb in combinations(0:(L-1), k)
        push!(ev, collect(comb))
    end
    ev
end
main()
