# =============================================================================
# imag_imgderive.jl —— img 映射的显式矩阵推导（N=1，8 维）
#
# ① 构造切片副本空间 ⊗_s C²（dim 2^k=8）上的 c_s/c†_s；
# ② img 弦规范候选（无弦 / 后缀 Z 弦 / 前缀 Z 弦 / (1−2n) 弦 × 顺序）对 17 探针标定；
# ③ 验证 img reproduces h-行（τ2 = Tr[img(S_K)·img(I)]）；
# ④ 解 ρ_I（64 维），检验 ρ_I ∈ img 像 span（= 非均匀失败的存在性判决）；
# ⑤ 若在 → 逆映射回多项式 = 闭式 I*。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_imgderive.jl
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

# ---------------- 几何（N=1 order-2 非均匀 [0.7,0.3]）----------------
const NL = 1
const CELLS = [0.7, 0.3]
const RLINKS = [0.7, 0.3]
const NL_INT = 2                     # 环区间数（nl；N=1 的 order-2 → 2 区间）
const LAT = GTEMPO.GrassmannLattice(N=NL_INT, δτ=BETA / NL_INT, contour=:imag, bands=1)
const L = length(LAT)                 # 6
const K = LAT.k                       # 3 切片
const DIM = 1 << K                    # 8
const EVENS = Vector{Vector{Int}}([Int[]])
for k in 2:2:L, comb in combinations(0:(L-1), k)
    push!(EVENS, collect(comb))
end
const NE = length(EVENS)
const CFGIDX = [1 + sum(1 << p for p in A; init=0) for A in EVENS]

# ---------------- 切片算符（8 维）----------------
c2 = ComplexF64[0 1; 0 0]; cd2 = ComplexF64[0 0; 1 0]; z2 = ComplexF64[1 0; 0 -1]
"""copy s（1..K）上的 2×2 嵌入。"""
function embed(m2::Matrix{ComplexF64}, s::Int)
    M = Matrix{ComplexF64}(I, DIM, DIM)
    for st in 0:(DIM-1)
        blk = m2[((st >> (s-1)) & 1) + 1, :]   # 行索引按该副本比特
        for out in 0:(DIM-1)
            (out >> (s-1)) & 1 == 0 || continue
            M[out + (1 << (s-1)) + 1, out + 1] = blk[2]
            M[out + 1, out + 1] = blk[1]
        end
    end
    return M
end
# 上面实现有误——改用标准 Kronecker 构造
function embed2(m2::Matrix{ComplexF64}, s::Int)
    factors = [s == j ? m2 : Matrix{ComplexF64}(I, 2, 2) for j in 1:K]
    M = factors[1]
    for j in 2:K
        M = kron(M, factors[j])
    end
    return M
end
const CS = [embed2(c2, s) for s in 1:K]; const CDS = [embed2(cd2, s) for s in 1:K]
const ZS = [embed2(z2, s) for s in 1:K]

# ---------------- img 弦规范候选 ----------------
# 生成器位置 pos（0-based L）：切片 s = 切片映射（pos p ↔ slice（K - p÷2）——由 lat 索引取）
slice_of(pos::Int) = K - (pos ÷ 2)            # pos 2m/2m+1 ↔ slice K-m
isbar_of(pos::Int) = (pos % 2) == 1           # 奇位 = bar

"""候选 img：单个生成器。string_mode: 0=无弦, 1=Z 后缀(>s), 2=Z 前缀(<s), 3=Z(≥s), 4=Z(≤s)"""
function img_gen(pos::Int, mode::Int)
    s = slice_of(pos)
    base = isbar_of(pos) ? CDS[s] : CS[s]
    if mode == 0
        return base
    elseif mode == 1
        W = s < K ? prod(ZS[s+1:K], init=Matrix{ComplexF64}(I, DIM, DIM)) : Matrix{ComplexF64}(I, DIM, DIM)
    elseif mode == 2
        W = s > 1 ? prod(ZS[1:s-1], init=Matrix{ComplexF64}(I, DIM, DIM)) : Matrix{ComplexF64}(I, DIM, DIM)
    elseif mode == 3
        W = prod(ZS[s:K], init=Matrix{ComplexF64}(I, DIM, DIM))
    else
        W = prod(ZS[1:s], init=Matrix{ComplexF64}(I, DIM, DIM))
    end
    return W * base
end

"""canonical 单项式（位集）→ img 算符（位置升序乘积）。"""
function img_mono(bits::Vector{Int}, mode::Int)
    M = Matrix{ComplexF64}(I, DIM, DIM)
    for p in sort(bits)
        M = M * img_gen(p, mode)
    end
    return M
end

# ---------------- 探针标定 ----------------
"""τ2 单链探针真值（来自 imag_trace_probe.jl，N=2 数据）——此处用 N=1 的
LAT 重新真测（注入 MPS + integrate）。"""
function probe_truth(pairsK::Vector{Tuple{Int,Int}}, pairsI::Vector{Tuple{Int,Int}})
    mk = GTEMPO.vacuumstate(ComplexF64, LAT)
    for (p1, p2) in sort(pairsK)
        mk = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mk)
    end
    mi = GTEMPO.vacuumstate(ComplexF64, LAT)
    for (p1, p2) in sort(pairsI)
        mi = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(1.0))), mi)
    end
    return real(GTEMPO.integrate(LAT, mk, mi))
end

function main()
    # 探针集（pos 0-based；L=6：0=s2? — 由 lat：slice s 的 eta pos = 2(K-s)，bar = +1）
    pos_eta(s) = 2 * (K - s); pos_bar(s) = pos_eta(s) + 1
    probes = [
        (Tuple{Int,Int}[], Tuple{Int,Int}[], 1.0),
        ([(pos_eta(3), pos_bar(3))], Tuple{Int,Int}[], 2.0),      # K=切片3 对 {0,1}
        ([(pos_eta(2), pos_bar(2))], Tuple{Int,Int}[], 2.0),
        ([(pos_eta(3), pos_bar(3))], [(pos_eta(3), pos_bar(3))], 3.0),
        ([(pos_eta(2), pos_bar(2))], [(pos_eta(3), pos_bar(3))], 4.0),
        ([(pos_eta(3), pos_bar(2))], Tuple{Int,Int}[], 1.0),      # 跨切片 {0,3}
        ([(pos_bar(3), pos_eta(2))], Tuple{Int,Int}[], 1.0),      # {1,2}
        ([(pos_eta(3), pos_bar(3)), (pos_eta(2), pos_bar(2))], Tuple{Int,Int}[], 4.0),
        ([(pos_eta(3), pos_bar(3)), (pos_eta(2), pos_bar(2))], [(pos_eta(1), pos_bar(1))], 8.0),
    ]
    println("img 规范标定（候选 × 探针）：")
    best = (-1.0, -1)
    for mode in 0:4
        err = 0.0
        perprobe = Float64[]
        for (kp, ip, want) in probes
            # 多项式先展开（θ²=0、canonical 符号），再线性映射 img
            M = Matrix{ComplexF64}(I, DIM, DIM)
            for (p1, p2) in kp
                monos = [(Int[], 1.0 + 0.0im), (sort([p1, p2]), 1.0 + 0.0im)]
                Mk = zeros(ComplexF64, DIM, DIM)
                for (loc, c) in monos
                    Mk .+= c .* img_mono(loc, mode)
                end
                M = M * Mk
            end
            Mi = Matrix{ComplexF64}(I, DIM, DIM)
            for (p1, p2) in ip
                monos = [(Int[], 1.0 + 0.0im), (sort([p1, p2]), 1.0 + 0.0im)]
                Mi2 = zeros(ComplexF64, DIM, DIM)
                for (loc, c) in monos
                    Mi2 .+= c .* img_mono(loc, mode)
                end
                Mi = Mi * Mi2
            end
            got = real(tr(M * Mi))
            err = max(err, abs(got - want))
            mode == 0 && push!(perprobe, got)
        end
        mode == 0 && println("  mode0 各探针 img 值：", perprobe)
        @printf("  mode %d : max err = %.3e\n", mode, err)
        err < 1e-10 && (best = (err, mode))
    end
    best[1] >= 0 && println("选中 mode = ", best[2])
    mode = best[2]

    # ---- img reproduces h-行？（τ2 = Tr[img(S_K)·img(I)]）
    println("\n验证：注入 I 单项式行的 img 重构")
    mk = GTEMPO.vacuumstate(ComplexF64, LAT)
    model = ToulouseIM(ϵ_d=0.7)
    sched = NativeExact.native_schedule(2, NL, BETA / NL)
    rls = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    for g in 1:NL
        s = g + 1
        pb = LAT[GrassmannContourIndex(s, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s, :τ, 1, false)]
        mk = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe;
                             coeff=ComplexF64(exp(-rls[g] * 0.7) - 1.0))), mk)
    end
    for p in 1:NL
        s1 = p + 1; s2 = p == NL ? 2 : p + 2
        pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
        pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
        mk = GTEMPO.apply!(exp(GTEMPO.GTerm(pb, pe; coeff=ComplexF64(-1.0))), mk)
    end
    p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
    p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
    mk = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mk)
    # img(S_K)：K 链多项式 → 单项式展开（用提取器思路：直接逐 GTerm 乘 img）
    rhoK = Matrix{ComplexF64}(I, DIM, DIM)
    for g in 1:NL
        s = g + 1
        pb = LAT[GrassmannContourIndex(s, :τ, 1, true)] - 1
        pe = LAT[GrassmannContourIndex(s, :τ, 1, false)] - 1
        c = exp(-rls[g] * 0.7) - 1.0
        rhoK = rhoK * (img_mono(sort([pe, pb]), mode) * c + Matrix{ComplexF64}(I, DIM, DIM))
    end
    for p in 1:NL
        s1 = p + 1; s2 = p == NL ? 2 : p + 2
        pb = LAT[GrassmannContourIndex(s1, :τ, 1, true)] - 1
        pe = LAT[GrassmannContourIndex(s2, :τ, 1, false)] - 1
        rhoK = rhoK * (img_mono(sort([pb, pe]), mode) * (-1.0) + Matrix{ComplexF64}(I, DIM, DIM))
    end
    rhoK = rhoK * (img_mono(sort([p1 - 1, p2 - 1]), mode) * (-1.0) + Matrix{ComplexF64}(I, DIM, DIM))
    # 对照：注入单项式 I 的真实 integrate vs Tr[rhoK·img(mono)]
    maxerr = 0.0
    for (ib, B) in enumerate(EVENS[1:min(20, NE)])
        mi = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (a, b2) in sort([(B[2i-1], B[2i]) for i in 1:(length(B) ÷ 2)])
            mi = GTEMPO.apply!(exp(GTEMPO.GTerm(a + 1, b2 + 1; coeff=ComplexF64(1.0))), mi)
        end
        ztrue = GTEMPO.integrate(LAT, mk, mi)
        Mi = length(B) == 0 ? Matrix{ComplexF64}(I, DIM, DIM) : img_mono(B, mode)
        zimg = real(tr(rhoK * Mi))
        maxerr = max(maxerr, abs(ztrue - zimg))
    end
    @printf("h-行 img 重构 max err = %.3e\n", maxerr)
end
main()
