# =============================================================================
# imag_o4split.jl —— order-4 残差的通道分离判决（h 行 vs g 行）
#
# 问题：联合行拟合残差 8e-5（非均匀 order-4）。分离拟合：
#   h-only（16 Z 行）残差 ~1e-12？ → I 槽对 Z 通道够用，坏在 g 通道（读出/插入约定）
#   h-only 也 8e-5？               → Z 通道本身结构性缺口（K 链或 I 槽表示）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4split.jl [tag...]
#   tag ∈ {n2o4ma, n1o4ma, n2o2ma}（默认三个全跑）
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

"""载入 h/g CSV；g CSV 每行 3 组 (z, ge, g1)。"""
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

function fit_residual(rowsV, rhs; ridge=1e-13)
    A = hcat(rowsV...)
    Ist = A * ((A' * A + ridge * I) \ rhs)
    return A' * Ist - rhs, Ist
end

function run_tag(outdir, tag, Nmacro, order)
    sched = NativeExact.native_schedule(order, Nmacro, BETA / Nmacro)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    L = 2 * (nl + 1)                       # imag lattice: k = N+1 节点 × (bar,eta)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    println("="^72)
    @printf("tag=%s  N=%d order=%d  nl=%d L=%d evensets=%d\n", tag, Nmacro, order, nl, L, NE)
    println("  cells = ", round.(sched.cells, digits=4))
    println("  ring  = ", round.(ring_links, digits=4))
    flush(stdout)

    hrows, grows = load_rows(outdir, tag, NE)
    heds = sort(collect(keys(hrows)))
    geds = sort(collect(keys(grows)))
    cells = sched.cells

    rowsH = [hrows[ed] for ed in heds]
    rhsH = [fock_anchor(cells, ring_links, ed)[1] for ed in heds]
    rowsG = [hrows[ed] .- grows[ed][NE .+ (1:NE)] .- real(fock_anchor(cells, ring_links, ed)[2]) .* hrows[ed] for ed in geds]
    rhsG = fill(0.0 + 0.0im, length(geds))

    rH, IstH = fit_residual(rowsH, rhsH)
    rG, IstG = fit_residual(rowsG, rhsG)
    rJ, IstJ = fit_residual(vcat(rowsH, rowsG), vcat(rhsH, rhsG))
    @printf("h-only 残差 |r| = %.3e   （16 Z 行）\n", norm(rH))
    @printf("g-only 残差 |r| = %.3e   （%d n 行）\n", norm(rG), length(geds))
    @printf("joint  残差 |r| = %.3e   （h+g 联合）\n", norm(rJ))

    println("\n[h-only 最大逐行残差]")
    perm = sortperm(abs.(rH), rev=true)
    for i in perm[1:min(4, end)]
        @printf("   Z@%-6.3f : %10.3e  (|rhs|=%.3e)\n", heds[i], abs(rH[i]), abs(rhsH[i]))
    end
    println("\n[g-only 最大逐行残差]")
    for (i, ed) in enumerate(geds)
        abs(rG[i]) > 1e-10 && @printf("   n@%-6.3f : %10.3e\n", ed, abs(rG[i]))
    end

    # 判决：用 h-only 解看 n 通道违反度（vend(I_h) vs (1−n)·hz(I_h)）
    println("\n[h-only 解的 n 通道检验]（Z 应精确、n 偏离=读出/插入问题）")
    maxv = 0.0
    for (i, ed) in enumerate(geds)
        hz = dot(hrows[ed], IstH)
        vend = dot(grows[ed][NE .+ (1:NE)], IstH)
        na = real(fock_anchor(cells, ring_links, ed)[2])
        npipe = real(1 - vend / hz)
        en = abs(npipe - na)
        ez = abs(real(hz) - real(rhsH[findfirst(==(ed), heds)]))
        maxv = max(maxv, en)
        @printf("   εd=%.3f : Zerr=%.2e  nerr=%.2e  (n=%.6f 锚 %.6f)\n", ed, ez, en, npipe, na)
    end
    @printf("h-only 解的最大 n 偏离 = %.3e\n", maxv)
    println(maxerr_summary(tag, norm(rH), norm(rG), norm(rJ), maxv))
    flush(stdout)
end

maxerr_summary(tag, rh, rg, rj, nv) =
    @sprintf("SUMMARY %s : h %.2e | g %.2e | joint %.2e | h-sol n %.2e",
             tag, rh, rg, rj, nv)

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    tags = isempty(ARGS) ? ["n2o2ma", "n1o4ma", "n2o4ma"] : ARGS
    specs = Dict("n2o2ma" => (2, 2), "n2o1ma" => (2, 1), "n3o2ma" => (3, 2),
                 "n1o4ma" => (1, 4), "n2o4ma" => (2, 4), "n1o2ma" => (1, 2))
    for tag in tags
        (N, o) = specs[tag]
        run_tag(outdir, tag, N, o)
    end
end
main()
