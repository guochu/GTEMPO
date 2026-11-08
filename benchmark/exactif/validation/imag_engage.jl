# =============================================================================
# imag_engage.jl —— Plan D 端到端落地：探针循环算法
#
# z(εd) = Σ_j c_j·integrate(mpsK(εd), probe_j)，探针 εd 无关（建一次复用）。
# n(εd) = 1 − vend/hz，vend 用末位插入通道 mpsK_end。
# 解 = imag_n2o4ma_dsol.csv（联合拟合输出）。对照 fock_anchor。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_engage.jl [εd...]
#   默认 εd = 0, 0.35, 0.7, 1.05, 1.4
# =============================================================================
using Printf
using LinearAlgebra
using Random
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation, TruncateDimCutoff

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

function lat_geom(nl::Int)
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    return LAT, L, EVENS
end

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

"""与 imag_probes.jl 完全一致的探针重建（同种子同顺序）。"""
function build_probes(LAT, EVENS, fam::String)
    rng = MersenneTwister(20261109)
    NE = length(EVENS)
    if fam == "s"
        probes = Vector{Tuple{Vector{Int}, Vector{Float64}}}()
        for B in EVENS
            np2 = length(B) ÷ 2
            σ = [rand(rng) < 0.5 ? -1.0 : 1.0 for _ in 1:np2]
            push!(probes, (B, σ))
        end
        return probes
    elseif fam == "p"
        probes = Vector{Tuple{Vector{Int}, Vector{Int}}}()
        for _ in 1:NE
            B1 = EVENS[rand(rng, 1:NE)]; B2 = EVENS[rand(rng, 1:NE)]
            push!(probes, (B1, B2))
        end
        return probes
    end
    error("fam")
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    NL, order = 2, 4
    sched = NativeExact.native_schedule(order, NL, BETA / NL)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    LAT, L, EVENS = lat_geom(nl)
    NE = length(EVENS)
    cells = sched.cells

    # 解文件
    sol = Vector{Tuple{String,Int,ComplexF64}}()
    for (i, ln) in enumerate(eachline(joinpath(outdir, "imag_n2o4ma_dsol.csv")))
        i == 1 && continue
        f = split(strip(ln), ",")
        push!(sol, (f[1], parse(Int, f[2]), ComplexF64(parse(Float64, f[3]), parse(Float64, f[4]))))
    end
    println("解项数：", length(sol)); flush(stdout)

    # 探针构建（按族缓存）
    inj = function(A::Vector{Int}, σ::Vector{Float64})
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (i, c) in zip(1:2:length(A), σ)
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(A[i] + 1, A[i+1] + 1; coeff=ComplexF64(c))), mps)
        end
        return mps
    end
    probecache = Dict{String, Any}()
    function getprobe(fam::String, jp::Int)
        cache = get!(probecache, fam) do
            Dict{Int, Any}()
        end
        get!(cache, jp) do
            if fam == "inj"
                A = EVENS[jp]
                inj(A, ones(Float64, length(A) ÷ 2))
            elseif fam == "s"
                probes = probecache["s_list"] === nothing ?
                    build_probes(LAT, EVENS, "s") : probecache["s_list"]
                (B, σ) = probes[jp]
                inj(B, σ)
            elseif fam == "p"
                probes = probecache["p_list"] === nothing ?
                    build_probes(LAT, EVENS, "p") : probecache["p_list"]
                (B1, B2) = probes[jp]
                m1 = inj(B1, ones(Float64, length(B1) ÷ 2))
                m2 = inj(B2, ones(Float64, length(B2) ÷ 2))
                GTEMPO.mult!(m1, m2; trunc=TruncateDimCutoff(256, 1e-13, 1))[1]
            end
        end
    end
    probecache["s_list"] = build_probes(LAT, EVENS, "s")
    probecache["p_list"] = build_probes(LAT, EVENS, "p")

    eds = isempty(ARGS) ? [0.0, 0.35, 0.7, 1.05, 1.4] :
        [parse(Float64, x) for x in ARGS]
    b1c = LAT[GrassmannContourIndex(1, conj=true, branch=:τ, band=1)]
    maxerrz = 0.0; maxerrn = 0.0
    for ed in eds
        mpsK = mk_kchain(LAT, ring_links, ed)
        tend = GTEMPO.GTerm(LAT[GrassmannContourIndex(LAT.Nτ, conj=false, branch=:τ, band=1)], b1c; coeff=1)
        mpsKend = GTEMPO._mult_A(tend, mpsK)
        z = 0.0 + 0.0im; vend = 0.0 + 0.0im
        for (k, (fam, jp, c)) in enumerate(sol)
            P = getprobe(fam, jp)
            z += c * GTEMPO.integrate(LAT, mpsK, P)
            vend += c * GTEMPO.integrate(LAT, mpsKend, P)
            k % 3000 == 0 && (println("  #$k/$(length(sol))"); flush(stdout))
        end
        npipe = real(1 - vend / z)
        (za, na) = fock_anchor(cells, ring_links, ed)
        ez = abs(real(z) - real(za)); en = abs(npipe - real(na))
        maxerrz = max(maxerrz, ez); maxerrn = max(maxerrn, en)
        @printf("εd=%-6.2f : Z = %12.8f（锚 %12.8f, err %.2e）  n = %9.6f（锚 %9.6f, err %.2e）\n",
                ed, real(z), real(za), ez, npipe, real(na), en)
        flush(stdout)
    end
    @printf("maxerr：Z %.2e / n %.2e  %s\n", maxerrz, maxerrn,
            (maxerrz < 5e-7 && maxerrn < 5e-7) ? "✅ PASS（order-4 探针循环精度）" : "✗")
end
main()
