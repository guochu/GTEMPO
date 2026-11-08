# =============================================================================
# imag_dfit.jl —— Plan D 判决：扩族拟合 + holdout 预测
#
# 列族：{inj(B)}（从 h 行 zeta 反演 raw）∪ {inj±(B,σ)}（s 族）∪ {mult!}（p 族）
# 行 = εd 锚点；rhs = fock_anchor Z。holdout = {0.35, 0.7, 1.05, 1.4}。
# 复合对 I-MPS 张量线性 → 任意解 c 可用 MPS 加法 mpsI = Σ c_j·probe_j 落地。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_dfit.jl
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

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    NL, order = 2, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links); L = 2 * (nl + 1)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    NE = length(EVENS)
    cells = sched.cells
    holdout = [0.35, 0.7, 1.05, 1.4]

    # ---- h 行（mobius 逆）→ zeta 正变换恢复 raw 注入行 ----
    hrows = Dict{Float64, Vector{ComplexF64}}()
    for file in sort(filter(f -> occursin("imag_n2o4ma_h", f) && endswith(f, ".csv"), readdir(outdir, join=true)))
        for (i, ln) in enumerate(eachline(file))
            i == 1 && continue
            f = split(strip(ln), ",")
            get!(hrows, parse(Float64, f[1]), zeros(ComplexF64, NE))[parse(Int, f[2])] =
                ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
        end
    end
    # raw_B = Σ_{sub ⊆ pairs(B)} h_{union(sub)}（zeta）
    rawrows = Dict{Float64, Vector{ComplexF64}}()
    for (ed, h) in hrows
        r = zeros(ComplexF64, NE)
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            np2 = length(B) ÷ 2
            for msk in Iterators.product(fill(0:1, np2)...)
                sub = Int[]
                for (i, bit) in enumerate(msk)
                    bit == 1 || continue
                    append!(sub, [B[2i-1], B[2i]])
                end
                val += h[findfirst(==(sort(sub)), EVENS)]
            end
            r[ib] = val
        end
        rawrows[ed] = r
    end
    println("raw 注入行恢复：$(length(rawrows)) εd"); flush(stdout)
    heds = sort(collect(keys(hrows)))
    holdout = [0.35, 0.7, 1.05, 1.4]
    anchors = setdiff(heds, holdout)
    println("锚点 $(length(anchors))（holdout $(length(holdout))）"); flush(stdout)

    # ---- 列装配（主 + holdout 同构） ----
    # 族数据载入：inj（rawrows）+ s/p CSV
    famdata = Dict{String, Dict{Float64, Vector{ComplexF64}}}()
    famnp = Dict{String, Int}()
    for fam in ("s", "p")
        data = Dict{Float64, Vector{ComplexF64}}()
        np = 0
        for part in (1, 2)
            file = joinpath(outdir, "imag_n2o4ma_$(fam)$(part).csv")
            isfile(file) || continue
            for (i, ln) in enumerate(eachline(file))
                i == 1 && continue
                f = split(strip(ln), ",")
                length(f) == 4 || (println("  族 ", fam, " part ", part, " 截断行 @", i, " 跳过"); break)
                ed = parse(Float64, f[1]); jp = parse(Int, f[2])
                g = get!(data, ed, zeros(ComplexF64, 0))
                if jp > length(g)
                    append!(g, zeros(ComplexF64, jp - length(g) + 1024)); data[ed] = g
                end
                g[jp] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
                np = max(np, jp)
            end
        end
        complete = np > 0 && all(ed -> haskey(data, ed), heds)
        if complete
            famdata[fam] = data; famnp[fam] = np
            println("族 ", fam, "：", np, " 列（完整）")
        else
            famnp[fam] = 0
            println("族 ", fam, "：数据不全，跳过")
        end
    end
    # 列布局：inj(NE) | s(np_s) | p(np_p)
    off_s = NE + 1
    off_p = NE + famnp["s"] + 1
    ncol = NE + famnp["s"] + famnp["p"]

    function probe_row(ed::Float64)
        r = Vector{ComplexF64}(undef, ncol)
        r[1:NE] .= rawrows[ed]
        for (fam, off) in (("s", off_s), ("p", off_p))
            np = famnp[fam]
            np == 0 && continue
            haskey(famdata[fam], ed) || error("族 $fam 缺 εd=$ed")
            r[off:off + np - 1] .= famdata[fam][ed][1:np]
        end
        return r
    end

    A = permutedims(hcat([probe_row(ed) for ed in anchors]...))
    rhs = [fock_anchor(cells, ring_links, ed)[1] for ed in anchors]
    @printf("设计矩阵 %d×%d\n", size(A, 1), size(A, 2)); flush(stdout)
    # 列归一化（压动态范围，不改变张成空间）
    cnrm = [norm(A[:, j]) for j in 1:size(A, 2)]
    cnrm[cnrm .== 0] .= 1.0
    An = A ./ cnrm'
    F = svd(An)
    println("σ 谱采样：", round.(F.S[[1, 2, 5, 10, 15, 20, 25, 30, 35, 40, end]], digits=6),
            "  rank(>1e-12·σ1) = ", count(>(1e-12 * F.S[1]), F.S))
    dtrue = norm(rhs - F.U * (F.U' * rhs))
    @printf("rhs 真距离 = %.3e（|rhs|=%.3e）\n", dtrue, norm(rhs)); flush(stdout)

    # ---- 联合拟合（IMAG_GFIT=1 时加载 g 通道行）----
    gfit = get(ENV, "IMAG_GFIT", "0") == "1"
    if gfit
        gdata = Dict{String, Dict{Float64, Vector{ComplexF64}}}()
        for fam in ("s", "p")
            famnp[fam] == 0 && continue
            data = Dict{Float64, Vector{ComplexF64}}()
            for part in (1, 2)
                file = joinpath(outdir, "imag_n2o4ma_$(fam)$(part)g.csv")
                isfile(file) || error("缺 $(file)")
                for (i, ln) in enumerate(eachline(file))
                    i == 1 && continue
                    f = split(strip(ln), ",")
                    length(f) == 4 || break
                    ed = parse(Float64, f[1]); jp = parse(Int, f[2])
                    g = get!(data, ed, zeros(ComplexF64, 0))
                    if jp > length(g)
                        append!(g, zeros(ComplexF64, jp - length(g) + 1024)); data[ed] = g
                    end
                    g[jp] = ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
                end
            end
            gdata[fam] = data
        end
        function gprobe_row(ed::Float64)
            r = Vector{ComplexF64}(undef, ncol)
            r[1:NE] .= rawrows[ed]   # 占位（inj 列在 g 行中不用——见下行覆盖）
            for (fam, off) in (("s", off_s), ("p", off_p))
                np = famnp[fam]
                np == 0 && continue
                haskey(gdata[fam], ed) || error("g 族 $fam 缺 εd=$ed")
                r[off:off + np - 1] .= gdata[fam][ed][1:np]
            end
            return r
        end
    end

    AH = permutedims(hcat([probe_row(ed) for ed in holdout]...))
    AHg = gfit ? permutedims(hcat([gprobe_row(ed) for ed in holdout]...)) : nothing
    println("\nkeep 扫描 + holdout 预测（Z）：")
    head = @sprintf("%-10s %10s", "keep", "res_anch")
    for ed in holdout; head *= @sprintf(" %11s", "Z@$ed"); end
    println(head)
    if gfit
        # 设计加倍：h 行（rhs=Z 锚）+ g 约束行（vend − (1−n)hz = 0）
        rows2 = Vector{Vector{ComplexF64}}(); rhs2 = ComplexF64[]
        for (i, ed) in enumerate(anchors)
            push!(rows2, A[i, :]); push!(rhs2, rhs[i])
            na = real(fock_anchor(cells, ring_links, ed)[2])
            prow = probe_row(ed); grow = gprobe_row(ed)
            push!(rows2, grow .- (1.0 - na) .* prow); push!(rhs2, 0.0 + 0.0im)
        end
        A2 = permutedims(hcat(rows2...))
        cnrm2 = [norm(A2[:, j]) for j in 1:size(A2, 2)]
        cnrm2[cnrm2 .== 0] .= 1.0
        An2 = A2 ./ cnrm2'
        global F, cnrm, An2, A2, rhs2
        F = svd(An2)
        cnrm = cnrm2
        A = A2; rhs = rhs2
        println("联合设计 ", size(A2, 1), "×", size(A2, 2),
                "  rank(>1e-12·σ1) = ", count(>(1e-12 * F.S[1]), F.S))
    end
    bestkeep = nothing; besterr = Inf
    for reltol in [1e-6, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14]
        keep = F.S .> reltol * F.S[1]
        Istn = F.V * ((F.U' * rhs) ./ ifelse.(keep, F.S, Inf))
        Ist = Istn ./ cnrm
        resa = norm(A * Ist - rhs)
        line = @sprintf("rel%-7.0e %10.2e", reltol, resa)
        herr = 0.0
        for (i, ed) in enumerate(holdout)
            pred = real(dot(AH[i, :] ./ cnrm, Istn))
            truth = real(fock_anchor(cells, ring_links, ed)[1])
            herr = max(herr, abs(pred - truth))
            line *= @sprintf(" %11.7f", pred)
        end
        if gfit
            for (i, ed) in enumerate(holdout)
                hz = dot(AH[i, :] ./ cnrm, Istn)
                vend = dot(AHg[i, :] ./ cnrm, Istn)
                npred = real(1 - vend / hz)
                truthn = real(fock_anchor(cells, ring_links, ed)[2])
                herr = max(herr, abs(npred - truthn))
                line *= @sprintf(" %9.6f", npred)
            end
        end
        line *= @sprintf("   maxhold=%.2e", herr)
        println(line); flush(stdout)
        if herr < besterr
            besterr = herr
            bestkeep = (reltol, Ist)
        end
    end
    println(gfit ? "真值Z                 2.5890924    2.2281719    1.9525846    1.7403380   （后 4 列 = n 预测，真值 0.419882/0.332466/0.249210/0.172405）"
                 : "真值                    2.5890924    2.2281719    1.9525846    1.7403380")
    # 解文件（best keep）
    (reltol, Ist) = bestkeep
    file = joinpath(outdir, "imag_n2o4ma_dsol.csv")
    open(file, "w") do io
        println(io, "fam,jp,re,im")
        for j in eachindex(Ist)
            abs(Ist[j]) > 1e-13 || continue
            fam = j <= NE ? "inj" : (j <= NE + famnp["s"] ? "s" : "p")
            jp = j <= NE ? j : (j <= NE + famnp["s"] ? j - NE : j - NE - famnp["s"])
            println(io, fam, ",", jp, ",", real(Ist[j]), ",", imag(Ist[j]))
        end
    end
    @printf("解已存 %s（best rel=%.0e, maxhold=%.2e，非零 %d 项）\n", file, reltol, besterr,
            count(x -> abs(x) > 1e-13, Ist))
    line0 = @sprintf("%-10s %10s", "真值", "")
    for ed in holdout
        line0 *= @sprintf(" %11.7f", real(fock_anchor(cells, ring_links, ed)[1]))
    end
    println(line0)
end
main()
