# =============================================================================
# imag_o4hold.jl —— Plan B 判决：锚点饱和的 holdout 预测收敛检验
#
# 原理：h 行 = canonical 基上的 h_b(K(εd))。函数族 F_resp 维度 R：当锚点数
#   > R 且采样在族上单射时，所有满足锚点的解的 h-函数相同 → min-norm 解的
#   网格外行为 = 真实行为。检验：拟合（排除 holdout 点）→ 预测 holdout 的
#   Z 值 → 随锚点增多/keep 扫描是否收敛到锚。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_o4hold.jl
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
    tag, NL, order = "n2o4ma", 2, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links); L = 2 * (nl + 1)
    EVENS = Vector{Vector{Int}}([Int[]])
    for m in 2:2:L, comb in combinations(0:(L-1), m)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    hrows, grows = load_rows(outdir, tag, NE)
    cells = sched.cells
    # g 行需要对应 εd 的 h 行（vend = (1−n)·hz 的 g 约束行）
    holdout = [0.35, 0.7, 1.05, 1.4]

    anchors = sort(collect(keys(hrows)))
    hanchors = setdiff(anchors, holdout)
    ganchors = sort(collect(keys(grows)))
    println("h 点：$(length(anchors))（holdout $(length(holdout))）  g 点：$(length(ganchors))")

    # 构造行（h 锚 + g 约束锚，holdout 全部排除）
    rowsV = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]; tagv = String[]
    for ed in hanchors
        push!(rowsV, hrows[ed]); push!(rhs, fock_anchor(cells, ring_links, ed)[1]); push!(tagv, "Z@$ed")
    end
    for ed in ganchors
        haskey(hrows, ed) || continue
        push!(rowsV, hrows[ed] .- grows[ed][NE .+ (1:NE)] .- real(fock_anchor(cells, ring_links, ed)[2]) .* hrows[ed])
        push!(rhs, 0.0 + 0.0im); push!(tagv, "n@$ed")
    end
    # holdout 行（只预测不拟合）
    rowsH = [hrows[ed] for ed in holdout]
    rhsH = [fock_anchor(cells, ring_links, ed)[1] for ed in holdout]

    A = hcat(rowsV...)
    AH = hcat(rowsH...)
    F = svd(A)
    @printf("锚点行矩阵 %d×%d，rank(σ>1e-12·σ₁)=%d，σ 谱末 6 个：%s\n",
            size(A, 2), size(A, 1), count(>(1e-12 * F.S[1]), F.S),
            string(round.(F.S[end-5:end], digits=12)))
    println("holdout 预测（真距离按 keep 截断扫描）：")
    header = @sprintf("%-10s %9s", "keep(tol)", "res_anch")
    for ed in holdout; header *= @sprintf(" %11s", "Z@$ed"); end
    println(header)
    for reltol in (1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12)
        keep = F.S .> reltol * F.S[1]
        Ist = F.U * ((F.V' * rhs) ./ ifelse.(keep, F.S, Inf))
        resa = norm(A' * Ist - rhs)
        line = @sprintf("%-10s %9.2e", "rel$reltol", resa)
        for (i, ed) in enumerate(holdout)
            pred = real(dot(rowsH[i], Ist))
            line *= @sprintf(" %11.7f", pred)
        end
        println(line)
    end
    line0 = @sprintf("%-10s %9s", "真值", "")
    for ed in holdout
        line0 *= @sprintf(" %11.7f", real(fock_anchor(cells, ring_links, ed)[1]))
    end
    println(line0)
    println("\n（若预测行随 keep 收敛到真值行 → 族饱和、min-norm 网格外正确 → 落地端到端）")
    # 额外：holdout 的 rank 分解——rhs_holdout 在列空间内的分量
    for (i, ed) in enumerate(holdout)
        vh = real(rhsH[i])
        best = Inf
        for reltol in (1e-7, 1e-9, 1e-11)
            keep = F.S .> reltol * F.S[1]
            Ist = F.U * ((F.V' * rhs) ./ ifelse.(keep, F.S, Inf))
            best = min(best, abs(real(dot(rowsH[i], Ist)) - vh))
        end
        @printf("holdout Z@%-5.2f = %.7f   最佳预测误差 = %.2e\n", ed, vh, best)
    end
end
main()
