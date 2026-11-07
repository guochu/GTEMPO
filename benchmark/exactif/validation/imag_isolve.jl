# =============================================================================
# imag_isolve.jl —— 从 h 网格精确解出所需 I 算符并检视其结构
#
# 方程：Σ_b I_b·h_εd(b) = Z_dense(εd)，16 个 εd × 16 个偶宇称未知数。
# I* 出来后与 row-exp 族（analytic_if）对照：
#   · 若 I* = exp-form → 读 Δ，拟合失败即 bug；
#   · 若 I* 非高斯 → 结构性结论：IF 需扩展（分辨率因子类），读出具体形式。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_isolve.jl
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
    Z = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z
    end
    return tr(Z)
end
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
    return -D
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
const NS = length(EVENS)
const IDXOF = Dict(A => i for (i, A) in enumerate(EVENS))

function load_h()
    hdata = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()
    for (i, ln) in enumerate(eachline(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_h.csv")))
        i == 1 && continue
        f = split(strip(ln), ",")
        key = (parse(Int, f[1]), parse(Float64, f[2]))
        ib = parse(Int, f[3])
        get!(hdata, key, zeros(ComplexF64, NS))[ib] = ComplexF64(parse(Float64, f[4]), parse(Float64, f[5]))
    end
    return hdata
end

function main()
    hdata = load_h()
    println("load_h keys = ", length(hdata))
    order = 2
    eds = sort(unique([k[2] for k in keys(hdata) if k[1] == order]))
    evenidx = [i for (i, A) in enumerate(EVENS) if true]   # 全部偶集合（每链偶宇称）
    println("εd 点数 = ", length(eds), "  偶集合 = ", NS)
    # 16×16：行 = εd，列 = 偶集合（即 I 多项式系数的支撑）
    n_use = min(length(eds), NS)
    Msys = [hdata[(order, eds[r])][c] for r in 1:n_use, c in 1:NS]
    zanch = [fock_anchor(CELLS, CELLS, eds[r]) for r in 1:n_use]
    # 解 I*（列维 NS，方程 n_use）：用全部 NS 列 + n_use 行；n_use=NS=32? NS=32
    println("方程数 n_use = $(n_use)，未知列 NS = $(NS)")
    # I 支撑限制到偶集合全集（32）——但多项式只取偶宇称集合（16）：odd 个腿的
    # 偶集合不存在——EVENS 全是偶腿集合 ✓ 32 个。方程 16 个 → 欠定；用
    # 高斯结构先验：先解 16×16（取前 16 个偶集合? 不行——需选支撑）。
    # 直接最小二乘/伪逆 + 基追踪风格：先试 pinv 最小范数解看结构。
    Istar = pinv(Msys) * zanch
    res = norm(Msys * Istar - zanch)
    @printf("pinv 最小范数解：残差 = %.3e\n", res)
    println("I* 非零项（偶集合 → 值）：")
    for (i, A) in enumerate(EVENS)
        abs(Istar[i]) > 1e-8 || continue
        @printf("   %-14s : %10.6f%s\n", A, real(Istar[i]),
                abs(imag(Istar[i])) > 1e-8 ? " + $(round(imag(Istar[i]), digits=6))im" : "")
    end
    println("（对照 row-exp 族：任意 Δ 的支撑 = ∅ + pairs(bar_i, eta_j)——6 个 pair 位）")
end
main()
