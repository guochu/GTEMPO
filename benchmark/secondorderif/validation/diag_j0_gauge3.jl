# ============================================================================
# diag_j0_gauge3.jl — τ 裸 gauge χ_bare 的通用公式推导（N=3，paritybits 基）
#
# 判决链（diag_j0_gauge.jl N=2）：补偿的作用 = 把 num 链 gauge 拉回 χ₀（无插入
# gauge）；χ_bare(topo) 用 21 特征（逐格带 paritybits + 全对乘积）GF(2) 拟合。
# 目标：找出 χ_required(topo) = χ₀ ⊕ χ_bare(topo) 的闭式 ⟹ 折算成 fold R 线
# (s0→tL) 逐项交叉补偿公式。
# 用法：julia --project=<repo> <repo>/benchmark/secondorderif/validation/diag_j0_gauge3.jl <repo>
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using GTEMPO
const TK = GTEMPO.TK
using Printf
using LinearAlgebra: mul!

ENV["TXCOMP"] = get(ENV, "GTX", "off")
const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
const N = 3

function coeffs(mps)
    L = length(mps)
    Ts = map(1:L) do i
        t = mps[i]
        d0i = TK.blockdim(TK.space(t, 1), TK.Z2Irrep(0)); d1i = TK.blockdim(TK.space(t, 1), TK.Z2Irrep(1))
        d0o = TK.blockdim(TK.domain(t)[1], TK.Z2Irrep(0)); d1o = TK.blockdim(TK.domain(t)[1], TK.Z2Irrep(1))
        T0 = zeros(ComplexF64, d0i + d1i, d0o + d1o); T1 = zeros(ComplexF64, d0i + d1i, d0o + d1o)
        for (f1, f2) in TK.fusiontrees(t)
            ql, p = Int(f1.uncoupled[1].n), Int(f1.uncoupled[2].n)
            qr = Int(f2.uncoupled[1].n)
            B = t[f1, f2]
            offq = ql == 0 ? 0 : d0i; offr = qr == 0 ? 0 : d0o
            for a in 1:size(B, 1), bb in 1:size(B, 3)
                v = B[a, 1, bb]
                iszero(v) && continue
                T = p == 0 ? T0 : T1
                T[offq + a, offr + bb] += v
            end
        end
        (T0, T1)
    end
    env = reshape(ComplexF64[1.0], 1, 1)
    for i in 1:L
        T0, T1 = Ts[i]
        n = size(env, 1)
        ne = zeros(ComplexF64, 2n, size(T1, 2))
        mul!(view(ne, 1:n, :), env, T0); mul!(view(ne, n+1:2n, :), env, T1)
        env = ne
    end
    co = GTEMPO.scaling(mps) .* env[:, 1]
    d = Dict{Int,ComplexF64}()
    for m in 0:(length(co) - 1)
        abs(co[m+1]) > 1e-10 && (d[m] = co[m+1])
    end
    return d
end

bath = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
mdl1 = ToulouseIM(ϵ_d = μ)
mdl2 = AndersonIM(ϵ_d = μ, U = 0.0)
c1 = WeylImagTime(Nτ = N, dτ = β / N, bands = 1)
c2 = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
dops1 = SecondOrderWeyl._fock_ops(1)
dops2 = SecondOrderWeyl._fock_ops(2)

function tauK(M::Int; ins::Vector{Tuple{Int,Matrix{ComplexF64}}})
    c = M == 1 ? c1 : c2
    wis = [WeylInsertion(:τ, b, mx) for (b, mx) in ins]
    system_factor(c, M == 1 ? mdl1 : mdl2; order = 2, coupling_norm = α, insertions = wis)
end

function sitemap(M::Int)
    lat = M == 1 ? weyl_lattice(c1) : weyl_lattice(c2)
    d = Dict{Int,Tuple{Int,Int,Bool}}()
    for j in 1:N, b in 1:M, cc in (false, true)
        s = lat[GrassmannContourIndex(j + 1, :τ, b, cc)]
        d[s - 1] = (j, b, cc)
    end
    return d
end
sm1 = sitemap(1); sm2 = sitemap(2)
# (cell, band) → [a 位, ā 位]（0-based）
siteab = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
for (s0, (j, b, cc)) in sm2
    a, ā = get(siteab, (j, b), (-1, -1))
    cc ? (siteab[(j, b)] = (a, s0)) : (siteab[(j, b)] = (s0, ā))
end
function embedmap(band::Int)
    em = Dict{Int,Int}()
    for (s2, (j, b, cc)) in sm2
        b == band || continue
        for (s1, (j1, _, cc1)) in sm1
            (j1 == j && cc1 == cc) && (em[s2] = s1; break)
        end
    end
    return em
end
em1 = embedmap(1); em2 = embedmap(2)
function splitmono(m::Int)
    ma = 0; mb = 0
    mm = m
    while mm != 0
        bt = mm & (-mm); mm -= bt
        s0 = trailing_zeros(bt)
        haskey(em1, s0) ? (ma |= 1 << em1[s0]) :
        haskey(em2, s0) ? (mb |= 1 << em2[s0]) : return nothing, nothing
    end
    return ma, mb
end

"cell j、band b 的 slot code paritybits（代码补偿的真实 c_b 基）"
slotpb(m, b, j) = begin
    a, ā = siteab[(j, b)]
    count_ones((((m >> a) & 1) ⊻ (((m >> ā) & 1) << 1))) % 2
end

d1 = Matrix{ComplexF64}(dops1[1]); dd1 = Matrix{ComplexF64}(dops1[1]')
D1 = Matrix{ComplexF64}(dops2[1]); DD1 = Matrix{ComplexF64}(dops2[1]')
D2 = Matrix{ComplexF64}(dops2[2]); DD2 = Matrix{ComplexF64}(dops2[2]')
n1ref = Matrix{ComplexF64}(dops1[1]' * dops1[1])
n1 = Matrix{ComplexF64}(dops2[1]' * dops2[1])

const FEATNAMES = vcat(["c$b@$k" for k in 1:N for b in 1:2],
                       [string("c", b, "@", k, "·c", b2, "@", k2)
                        for k in 1:N for b in 1:2 for k2 in (k+1):N for b2 in 1:2])
nfeat = length(FEATNAMES)
function feats(m)
    cb = [slotpb(m, b, k) for k in 1:N for b in 1:2]
    idx(k, b) = (k - 1) * 2 + b
    cross = [cb[idx(k, b)] & cb[idx(k2, b2)] for k in 1:N for b in 1:2
             for k2 in (k+1):N for b2 in 1:2]
    return vcat(cb, cross)
end

"gauge 测量 + GF(2) 拟合，返回 χ 的特征系数 λ"
function gauge(topo::String, insA, insB, ins2; show::Bool = true)
    Pa = coeffs(tauK(1; ins = insA))
    Pb = coeffs(tauK(1; ins = insB))
    Q = coeffs(tauK(2; ins = ins2))
    signs = Dict{Int,Int}()
    ncmp = 0; nstruct = 0
    structex = Tuple{Int,Int,ComplexF64,ComplexF64,ComplexF64}[]
    for (m, vq) in Q
        ma, mb = splitmono(m)
        ma === nothing && continue
        va = get(Pa, ma, 0.0im); vb = get(Pb, mb, 0.0im)
        (abs(va) > 1e-9 && abs(vb) > 1e-9 && abs(vq) > 1e-12) || continue
        r = vq / (va * vb)
        ncmp += 1
        if abs(abs(r) - 1) < 1e-6
            signs[m] = real(r) > 0 ? 1 : -1
        else
            nstruct += 1
            nstruct <= 8 && push!(structex, (m, mb, vq, va * vb, r))
        end
    end
    isempty(structex) || begin
        println("  结构样本（mask, mb(带2参照位), Q, Pa·Pb, ρ）：")
        for (m, mb, vq, vp, r) in structex
            @printf("    mask %s | mb %s | Q=%+.5f%+.5fi | P=%+.5f%+.5fi | ρ=%+.4f%+.4fi |q|/|p|=%.4f\n",
                    string(m, base=2, pad=12), string(mb, base=2, pad=6),
                    real(vq), imag(vq), real(vp), imag(vp), real(r), imag(r), abs(vq)/abs(vp))
        end
    end
    ms = collect(keys(signs))
    rows = [feats(m) for m in ms]
    b = [signs[m] == -1 ? 1 : 0 for m in ms]
    piv = Tuple{Int,Int}[]; r0 = 1
    for c in 1:nfeat
        rp = findfirst(r -> r >= r0 && rows[r][c] == 1, r0:length(ms))
        rp === nothing && continue
        rp += r0 - 1
        rows[rp], rows[r0] = rows[r0], rows[rp]
        b[rp], b[r0] = b[r0], b[rp]
        for r in 1:length(rows)
            if r != r0 && rows[r][c] == 1
                rows[r] .⊻= rows[r0]
                b[r] ⊻= b[r0]
            end
        end
        push!(piv, (r0, c)); r0 += 1
        r0 > length(rows) && break
    end
    incons = count(r -> all(iszero, rows[r]) && b[r] != 0, 1:length(rows))
    λ = zeros(Int, nfeat)
    for (r, c) in piv
        λ[c] = b[r]
    end
    show && begin
        act = [FEATNAMES[x] for x in 1:nfeat if λ[x] == 1]
        @printf("[%s] 可比 %d；|ρ|≠1 = %d；-1 = %d；拟合（不一致行 %d）： %s\n",
                topo, ncmp, nstruct, count(==(-1), values(signs)), incons,
                isempty(act) ? "（恒 +1）" : join(act, " ⊕ "))
    end
    return λ, incons
end

println("== N=3 TXCOMP=$(ENV["TXCOMP"]) gauge（paritybits 基）==")
const NI = Tuple{Int,Matrix{ComplexF64}}[]
λ0, _ = gauge("无插入", NI, NI, NI)
topos = Tuple{String,Vector{Tuple{Int,Matrix{ComplexF64}}}}[
    ("(1,0)", [(1, D1), (0, DD1)]),
    ("(2,0)", [(2, D1), (0, DD1)]),
    ("(3,0)", [(3, D1), (0, DD1)]),
    ("(2,1)", [(2, D1), (1, DD1)]),
    ("(3,1)", [(3, D1), (1, DD1)]),
    ("(3,2)", [(3, D1), (2, DD1)]),
    ("(1,1)", [(1, D1), (1, DD1)]),
    ("(3,3)", [(3, D1), (3, DD1)]),
    ("occ(0)", [(0, n1)]),
]
results = Dict{String,Tuple{Vector{Int},Int}}()
for (lbl, ins) in topos
    results[lbl] = gauge(lbl, NI, NI, ins)
end
println("\n== χ_required = χ₀ ⊕ χ_bare（每拓扑的特征差分）==")
for (lbl, _) in topos
    λb, ic = results[lbl]
    ic > 0 && (println("  $lbl: 拟合不一致（$ic 行），跳过"); continue)
    d = λ0 .⊻ λb
    act = [FEATNAMES[x] for x in 1:nfeat if d[x] == 1]
    println("  $lbl: ", isempty(act) ? "（= χ₀，无需补偿）" : join(act, " ⊕ "))
end
println("diag_j0_gauge3 完成")
