# =============================================================================
# imag_bilin.jl —— 位置模代数/Λ 泛函模型在 order-4 上的判决（nl=6, L=14）
#
# 模型：级联积 = Grassmann 多项式积；复合 integrate(K,I) = Λ(P_K × P_I)。
#   Λ = 线性泛函（偶单项式上值，真空行 Möbius 读出）。
# 判决：
#   V1：apply! 级联 == 多项式积？（mps_to_dense(mpsK) vs poly_mul 因子积）
#   V2：raw 行 == Λ(P_K × P_I(B))？（全 8192 B 对比）
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_bilin.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation, Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

"""θ_a·θ_b 的 shuffle 符号。"""
function gr_sign(a::Vector{Int}, b::Vector{Int})
    cnt = 0
    for p in a, q in b
        p > q && (cnt += 1)
    end
    return (-1.0)^cnt
end

"""多项式 Grassmann 积（bitmask+1 索引的向量）。"""
function poly_mul(u::Vector{ComplexF64}, v::Vector{ComplexF64}, L::Int)
    out = zeros(ComplexF64, 1 << L)
    for ia in eachindex(u)
        ua = u[ia]; ua == 0 && continue
        (ia - 1) & 0 == 0 || true
        for ib in eachindex(v)
            vb = v[ib]; vb == 0 && continue
            (ia - 1) & (ib - 1) != 0 && continue
            A = Int[p for p in 0:(L-1) if (ia - 1) >> p & 1 == 1]
            B = Int[p for p in 0:(L-1) if (ib - 1) >> p & 1 == 1]
            out[((ia - 1) | (ib - 1)) + 1] += ua * vb * gr_sign(A, B)
        end
    end
    return out
end

function lat_geom(nl::Int)
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    return LAT, L, EVENS
end

"""pipeline-v2 K 链（imag_general 同款）。"""
function mk_kchain(LAT, ring_links0::Vector{Float64}, eps_d::Float64)
    nl = length(ring_links0)
    mps = GTEMPO.vacuumstate(ComplexF64, LAT)
    for p in 1:nl
        iszero(eps_d) || begin
            sg = p + 1
            gb = LAT[GrassmannContourIndex(sg, :τ, 1, true)]
            ge = LAT[GrassmannContourIndex(sg, :τ, 1, false)]
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(gb, ge;
                                 coeff=ComplexF64(exp(-ring_links0[p] * eps_d) - 1.0))), mps)
        end
        s1 = p + 1
        s2 = p == nl ? 2 : p + 2
        rb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
        re = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(rb, re; coeff=ComplexF64(-1.0))), mps)
    end
    p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
    p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
    mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mps)
    return mps
end

"""GMPS → canonical 稠密系数表。"""
function mps_to_dense(x::GTEMPO.GrassmannMPS, L::Int)
    n = length(x)
    n == L || throw(ArgumentError("site 数应等于 L，n=$n L=$L"))
    A = ones(ComplexF64, 1, 1)
    for j in 1:n
        t = x[j]
        spL, domR = codomain(t), domain(t)
        dL, dR = spL[1], domR[1]
        dimL = [blockdim(dL, Z2Irrep(q)) for q in 0:1]
        dimR = [blockdim(dR, Z2Irrep(q)) for q in 0:1]
        offR = (0, dimR[1])
        newA = zeros(ComplexF64, dimR[1] + dimR[2], 2 * size(A, 2))
        for (f1, f2) in fusiontrees(t)
            qL = f1.uncoupled[1].n; bp = f1.uncoupled[2].n; qR = f2.uncoupled[1].n
            m = t[f1, f2]
            for α in 1:dimL[qL+1], β in 1:dimR[qR+1], cfg in 1:size(A, 2)
                newA[offR[qR+1] + β, cfg + (bp << (j - 1))] +=
                     m[α, 1, β] * A[α + (qL == 1 ? dimL[1] : 0), cfg]
            end
        end
        A = newA
    end
    return vec((scaling(x)^L) .* A)
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    NL, order = 2, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    LAT, L, EVENS = lat_geom(nl)
    NE = length(EVENS)
    EVIDX = fill(0, 1 << L)
    for (ib, A) in enumerate(EVENS)
        EVIDX[1 + sum(1 << p for p in A; init=0)] = ib
    end
    idxof(A) = EVIDX[1 + sum(1 << p for p in A; init=0)]
    println("nl=$nl L=$L evensets=$NE dim=2^L=$((1 << L))"); flush(stdout)

    vac = GTEMPO.vacuumstate(ComplexF64, LAT)
    injected = function(A::Vector{Int})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for i in 1:2:length(A)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(A[i] + 1, A[i+1] + 1; coeff=ComplexF64(1.0))), mps)
        end
        return mps
    end

    # ---------------- 真空行 → Λ ----------------
    Λfile = joinpath(outdir, "imag_n2o4ma_lam.csv")
    if !isfile(Λfile)
        println("测量真空行..."); flush(stdout)
        raw = Vector{ComplexF64}(undef, NE)
        t0m = time()
        for (ib, A) in enumerate(EVENS)
            ts = time()
            mpsB = injected(A)
            t1 = time()
            raw[ib] = GTEMPO.integrate(LAT, vac, mpsB)
            t2 = time()
            ib <= 3 && (println("  B#$(ib): inj $(round(t1-ts, digits=2))s integ $(round(t2-t1, digits=2))s"); flush(stdout))
            ib % 500 == 0 && (println("  ...$ib/$NE ($(round(time()-t0m, digits=1))s)"); flush(stdout))
        end
        Λ = zeros(ComplexF64, NE)
        for (ib, B) in enumerate(EVENS)
            val = 0.0 + 0.0im
            np2 = length(B) ÷ 2
            for msk in Iterators.product(fill(0:1, np2)...)
                sub = Int[]
                for (i, bit) in enumerate(msk)
                    bit == 1 || continue
                    append!(sub, [B[2i-1], B[2i]])
                end
                sgn = (-1.0)^((length(B) - length(sub)) ÷ 2)
                val += sgn * raw[idxof(sort(sub))]
            end
            Λ[ib] = val
        end
        open(Λfile, "w") do io
            println(io, "ib,re,im")
            for (ib, A) in enumerate(EVENS)
                abs(Λ[ib]) > 1e-12 || continue
                println(io, ib, ",", real(Λ[ib]), ",", imag(Λ[ib]))
            end
        end
        println("Λ 已存（非零 ", count(x -> abs(x) > 1e-12, Λ), "/$NE）")
    end
    Λ = zeros(ComplexF64, NE)
    for (i, ln) in enumerate(eachline(Λfile))
        i == 1 && continue
        f = split(strip(ln), ",")
        Λ[parse(Int, f[1])] = ComplexF64(parse(Float64, f[2]), parse(Float64, f[3]))
    end
    println("Λ 载入完成（非零 ", count(x -> abs(x) > 1e-12, Λ), "）"); flush(stdout)

    # Λ → 全多项式向量（bitmask+1 索引）
    Λfull = zeros(ComplexF64, 1 << L)
    for (ib, A) in enumerate(EVENS)
        abs(Λ[ib]) > 1e-12 || continue
        Λfull[1 + sum(1 << p for p in A; init=0)] = Λ[ib]
    end

    # ---------------- V1：apply! == 多项式积？ ----------------
    println("\n== V1：K 链提取 vs 因子多项式积 ==")
    for ed in [0.0, 0.7]
        mpsK = mk_kchain(LAT, ring_links, ed)
        S_mps = mps_to_dense(mpsK, L)
        # 因子积（同 mk_kchain 顺序）
        P = zeros(ComplexF64, 1 << L); P[1] = 1.0
        for p in 1:nl
            iszero(ed) || begin
                sg = p + 1
                gb = LAT[GrassmannContourIndex(sg, :τ, 1, true)] - 1
                ge = LAT[GrassmannContourIndex(sg, :τ, 1, false)] - 1
                F = zeros(ComplexF64, 1 << L); F[1] = 1.0
                F[1 + (1 << gb) | (1 << ge)] = exp(-ring_links[p] * ed) - 1.0
                P = poly_mul(F, P, L)
            end
            s1 = p + 1
            s2 = p == nl ? 2 : p + 2
            rb = LAT[GrassmannContourIndex(s1, :τ, 1, true)] - 1
            re = LAT[GrassmannContourIndex(s2, :τ, 1, false)] - 1
            F = zeros(ComplexF64, 1 << L); F[1] = 1.0
            F[1 + (1 << rb) | (1 << re)] = -1.0
            P = poly_mul(F, P, L)
        end
        p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)] - 1
        p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)] - 1
        F = zeros(ComplexF64, 1 << L); F[1] = 1.0
        F[1 + (1 << p1) | (1 << p2)] = -1.0
        P = poly_mul(F, P, L)
        err = maximum(abs.(S_mps - P))
        @printf("  εd=%.2f : max|S_mps − P_级联积| = %.3e   %s\n", ed, err,
                err < 1e-10 ? "✅ apply!=多项式积" : "✗ 表示差异")
        println("  V1 错配单项式（b：S_mps vs P_级联积）：")
        d1 = [(abs(S_mps[i] - P[i]), i) for i in eachindex(P)]
        sort!(d1, by = x -> -x[1])
        for (d, i) in d1[1:min(8, count(x -> x[1] > 1e-9, d1))]
            bits = [p for p in 0:(L-1) if (i - 1) >> p & 1 == 1]
            @printf("    b=%-22s : S=%9.4f  P=%9.4f\n", isempty(bits) ? "∅" : join(bits, ","), real(S_mps[i]), real(P[i]))
        end
    end

    # ---------------- V2：h 行 == G(b) = Σ_a PK_a·gr_sign(a,b)·Λ_{a∪b} ----------------
    # raw_B = Σ_{sub⊆pairs(B)} G(sub)（级联展开系数 +1），mobius 逆给出 h_row = G(b)。
    # 故直接对照 45 个已有 CSV 行，零新测量。
    println("\n== V2：h 行 vs G(b) 预测 =="); flush(stdout)
    L14 = 1 << L
    bits_of = Vector{Vector{Int}}(undef, L14)
    for m in 0:(L14 - 1)
        bits_of[m + 1] = Int[p for p in 0:(L-1) if m >> p & 1 == 1]
    end
    Λidx = [i - 1 for i in eachindex(Λfull) if abs(Λfull[i]) > 1e-14]
    println("Λ 支集大小：", length(Λidx))
    println("Λ 支集：", length(Λidx)); flush(stdout)

    hrows = Dict{Float64, Vector{ComplexF64}}()
    for file in sort(filter(f -> occursin("imag_n2o4ma_h", f) && endswith(f, ".csv"), readdir(outdir, join=true)))
        for (i, ln) in enumerate(eachline(file))
            i == 1 && continue
            f = split(strip(ln), ",")
            get!(hrows, parse(Float64, f[1]), zeros(ComplexF64, NE))[parse(Int, f[2])] =
                ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))
        end
    end
    heds = sort(collect(keys(hrows)))
    println("对照 h 行：$(length(heds)) 个 εd"); flush(stdout)

    for ed in heds
        mpsK = mk_kchain(LAT, ring_links, ed)
        PK = mps_to_dense(mpsK, L)
        G = zeros(ComplexF64, NE)
        for (ib, B) in enumerate(EVENS)
            mb = sum(1 << p for p in B; init=0)
            Bv = bits_of[mb + 1]
            acc = 0.0 + 0.0im
            for li in Λidx
                (li & mb) != mb && continue      # λ ⊇ b
                a = li & ~mb                      # a = λ \ b
                av = bits_of[a + 1]
                cnt = 0
                for pp in av, qq in Bv
                    pp > qq && (cnt += 1)
                end
                acc += Λfull[li + 1] * PK[a + 1] * (cnt % 2 == 0 ? 1.0 : -1.0)
            end
            G[ib] = acc
        end
        err = maximum(abs.(G - hrows[ed]))
        @printf("  εd=%-6.3f : max|G − h_row| = %.3e %s\n", ed, err,
                err < 1e-9 ? "✅" : "✗")
        ed == heds[1] || (flush(stdout); continue)
        println("  最大错配单项式（b：G 预测 vs h 实测）：")
        diffs = [(abs(G[i] - hrows[ed][i]), i) for i in 1:NE]
        sort!(diffs, by = x -> -x[1])
        for (d, i) in diffs[1:min(10, count(x -> x[1] > 1e-9, diffs))]
            @printf("    b=%-22s : G=%9.4f  h=%9.4f\n", join(EVENS[i], ","), real(G[i]), real(hrows[ed][i]))
        end
        flush(stdout)
    end
end
main()
