# =============================================================================
# imag_istar.jl —— 全系数 I** 精确解与结构读取
#
# 行：Z-行（h，rank 3）+ n-行（变体 1−g_end − n·h，rank ~8）。32 个偶宇称
# 未知数，min-norm 解。输出全部非零项与"pair-exp 预测值"的偏差表。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_istar.jl
# =============================================================================
using Printf
using LinearAlgebra
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
const L = length(LAT)
const EVENS = Vector{Vector{Int}}([Int[]])
for k in (2, 4, 6), comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end

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
    # 用全部 h-εd 行（rank 3）+ 全部 n-行
    hed_all = sort([k[2] for k in keys(hdata) if k[1] == 2])
    rows = Vector{ComplexF64}[]; rhs = ComplexF64[]
    for ed in hed_all
        push!(rows, hdata[(2, ed)][:, 1]); push!(rhs, fock_anchor(CELLS, CELLS, ed)[1])
    end
    for (r, ed) in enumerate(eds)
        hz = hdata[(2, ed)][:, 1]
        gv = gdata[(2, ed)]
        vrow = hz .- gv[:, 2]                       # 1 − g_end 变体
        nanc = real(fock_anchor(CELLS, CELLS, ed)[2])
        push!(rows, vrow .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
    end
    A = hcat(rows...)                                # 32 × nrows（列 = 行向量）
    b = rhs
    # 约束：Aᵀ I = b（nrows 行），min-norm：I = A · (AᵀA)⁻¹ b —— A: 32×nrows
    Ist = A * ((A' * A) \ b)
    res = norm(A' * Ist - b)
    @printf("约束行数 = %d，残差 = %.3e\n", length(rhs), res)
    println("\nI** 全部非零项（|v|>1e-8）：")
    for (i, Aset) in enumerate(EVENS)
        v = Ist[i]
        abs(v) > 1e-8 || continue
        @printf("   %-16s : %11.7f%s\n", Aset, real(v),
                abs(imag(v)) > 1e-8 ? " + $(round(imag(v), digits=7))im" : "")
    end
    println("\n行空间奇异性检查：svd(A' A) 最小奇异值 = ", minimum(svd(A' * A).S |> collect))
end
main()
