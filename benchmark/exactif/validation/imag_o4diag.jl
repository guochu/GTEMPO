# =============================================================================
# imag_o4diag.jl —— order-4 联合行不一致的残差解剖
#
# 输入：imag_n2o4ma_{h,g}{1,2}.csv（nl 修复后的行数据）
# 分析：
#   1. A 的 SVD 谱与有效秩；rhs 到 range(A) 的距离分解
#   2. min-norm 解的逐行残差（Z-行 vs n-行 × εd）
#   3. 残差向量与 h/g 行结构的关联（哪个通道、哪个 εd 段）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4diag.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const NL = 2
const WB = 1.0
const AB = 0.5
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

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    # 几何（order 4: nl = 3N = 6, L = 12）
    sched = NativeExact.native_schedule(4, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    @printf("nl=%d L=%d evensets=%d\n", nl, L, NE)

    hrows = Dict{Float64, Vector{ComplexF64}}()
    for part in (1, 2)
        for (i, ln) in enumerate(eachline(joinpath(outdir, "imag_n2o4ma_h$part.csv")))
            i == 1 && continue
            f = split(strip(ln), ",")
            get!(hrows, parse(Float64, f[1]), zeros(ComplexF64, NE))[parse(Int, f[2])] =
                ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
        end
    end
    grows = Dict{Float64, Vector{ComplexF64}}()
    for part in (1, 2)
        for (i, ln) in enumerate(eachline(joinpath(outdir, "imag_n2o4ma_g$part.csv")))
            i == 1 && continue
            f = split(strip(ln), ",")
            g = get!(grows, parse(Float64, f[1]), zeros(ComplexF64, 3NE))
            ib = parse(Int, f[2])
            g[ib] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
            g[NE + ib] = ComplexF64(parse(Float64, f[5]), parse(Float64, f[6]))
            g[2NE + ib] = ComplexF64(parse(Float64, f[7]), parse(Float64, f[8]))
        end
    end

    heds = sort(collect(keys(hrows)))
    geds = sort(collect(keys(grows)))
    cells = sched.cells
    rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]; rowtag = String[]
    for ed in heds
        push!(rows, hrows[ed]); push!(rhs, fock_anchor(cells, ring_links, ed)[1])
        push!(rowtag, "Z@$(ed)")
    end
    for ed in geds
        hz = hrows[ed]
        vend = grows[ed][NE .+ (1:NE)]
        nanc = real(fock_anchor(cells, ring_links, ed)[2])
        push!(rows, hz .- vend .- nanc .* hz); push!(rhs, 0.0 + 0.0im)
        push!(rowtag, "n@$(ed)")
    end
    A = hcat(rows...)                       # NE × nrows
    @printf("A: %d×%d\n", size(A, 1), size(A, 2))

    # SVD 谱
    F = svd(A)
    sv = F.S
    @printf("A 奇异值（前 12）：")
    println(round.(sv[1:min(12, length(sv))], digits=6))
    for tol in (1e-14, 1e-10, 1e-6)
        @printf("  rank(tol=%.0e) = %d\n", tol, count(>(tol * sv[1]), sv))
    end

    # min-norm 解与逐行残差
    Ist = A * ((A' * A + 1e-12 * I) \ rhs)
    r = A' * Ist - rhs
    @printf("\n联合残差 |r| = %.3e\n", norm(r))
    println("逐行残差：")
    for (i, tg) in enumerate(rowtag)
        abs(r[i]) > 1e-9 && @printf("   %-10s : %10.3e  (|rhs|=%.3e)\n", tg, abs(r[i]), abs(rhs[i]))
    end

    # rhs 到 range(A) 的距离（残差即该距离——分解到哪些奇异方向）
    proj = F.U * (F.U' * rhs)               # rhs 在列空间（NE 维）上的投影——注意 A 列空间
    # 直接：残差 r 是否属于 A 的"左零空间"分量
    # r = rhs − A'A(A'A)⁻¹rhs 的一部分……直接看 r 与各左奇异向量的内积
    Av = svd(A' * A)
    println("\n残差分解（|⟨u_k, rhs⟩| 前 10，u_k = A'A 的左奇异向量）：")
    for k in 1:min(10, length(Av.S))
        val = abs(dot(Av.U[:, k], rhs))
        val > 1e-8 && @printf("   σ = %10.3e : |⟨u,rhs⟩| = %10.3e\n", Av.S[k], val)
    end
end
main()
