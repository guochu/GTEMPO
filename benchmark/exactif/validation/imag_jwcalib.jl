# =============================================================================
# imag_jwcalib.jl —— JW 弦 img 同态标定
#
# img(θ_p) = (JW弦)·(c/c†)_{slice(p)}，弦方向/范围 4 候选 + 无弦对照。
# 标定：①17 探针（N=2 lattice 真实 integrate）；②同态性 img(A·B) == img(A)img(B)；
# ③锚点 ⟨0|img(K(0))|0⟩ == 真实 integrate(mpsK(0), vacuum)。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_jwcalib.jl
# =============================================================================
using Printf
using LinearAlgebra
using Combinatorics: combinations
using GTEMPO
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const NL = 2
const WB, AB = 1.0, 0.5
HB = ComplexF64[WB;;]
VB = ComplexF64[sqrt(AB) + 0.0im]
hcoupled() = NativeExact._coupled_h(VB, HB)

const LAT = GTEMPO.GrassmannLattice(N=NL, δτ=BETA / NL, contour=:imag, bands=1)
const L = length(LAT)
const K = LAT.k
const DIM = 1 << K
const VAC = begin v = zeros(ComplexF64, DIM); v[1] = 1.0; v end

c2 = ComplexF64[0 1; 0 0]; cd2 = ComplexF64[0 0; 1 0]; z2 = ComplexF64[1 0; 0 -1]
function embed2(m2::Matrix{ComplexF64}, s::Int)
    factors = [s == j ? m2 : Matrix{ComplexF64}(I, 2, 2) for j in 1:K]
    M = factors[1]
    for j in 2:K
        M = kron(M, factors[j])
    end
    return M
end
const CS = [embed2(c2, s) for s in 1:K]
const CDS = [embed2(cd2, s) for s in 1:K]
const ZS = [embed2(z2, s) for s in 1:K]

slice_of(pos::Int) = K - (pos ÷ 2)
isbar_of(pos::Int) = (pos % 2) == 1

"""JW 弦 img。mode: 0 无弦；1 Z(>s)；2 Z(<s)；3 Z(≥s)；4 Z(≤s)；5 Z(≠s)。"""
function img_gen(pos::Int, mode::Int)
    s = slice_of(pos)
    base = isbar_of(pos) ? CDS[s] : CS[s]
    I8 = Matrix{ComplexF64}(I, DIM, DIM)
    if mode == 0
        return base
    elseif mode == 1
        return s < K ? prod(ZS[s+1:K]; init=I8) * base : base
    elseif mode == 2
        return s > 1 ? prod(ZS[1:s-1]; init=I8) * base : base
    elseif mode == 3
        return prod(ZS[s:K]; init=I8) * base
    elseif mode == 4
        return prod(ZS[1:s]; init=I8) * base
    else
        W = Matrix{ComplexF64}(I, DIM, DIM)
        for t in 1:K
            t == s || (W *= ZS[t])
        end
        return W * base
    end
end
img_mono(bits::Vector{Int}, mode::Int) = begin
    M = Matrix{ComplexF64}(I, DIM, DIM)
    for p in sort(bits)
        M = M * img_gen(p, mode)
    end
    return M
end
vevc(M) = ComplexF64(dot(VAC, M * VAC))

function main()
    # ---------- ① 探针标定（17 探针真值用真实 integrate 现测）----------
    function inj(pairs::Vector{Tuple{Int,Int}})
        m = GTEMPO.vacuumstate(ComplexF64, LAT)
        for (a, b) in sort(pairs)
            m = GTEMPO.apply!(exp(GTEMPO.GTerm(a + 1, b + 1; coeff=ComplexF64(1.0))), m)
        end
        return m
    end
    vac = GTEMPO.vacuumstate(ComplexF64, LAT)
    pe(s) = 2 * (K - s); pb(s) = pe(s) + 1
    probes = [
        (Tuple{Int,Int}[], Tuple{Int,Int}[]),
        ([(pe(3), pb(3))], Tuple{Int,Int}[]),
        ([(pe(2), pb(2))], Tuple{Int,Int}[]),
        ([(pe(3), pb(3))], [(pe(3), pb(3))]),
        ([(pe(2), pb(2))], [(pe(3), pb(3))]),
        ([(pe(3), pb(2))], Tuple{Int,Int}[]),
        ([(pb(3), pe(2))], Tuple{Int,Int}[]),
        ([(pe(3), pb(3)), (pe(2), pb(2))], Tuple{Int,Int}[]),
        ([(pe(3), pb(3)), (pe(2), pb(2))], [(pe(1), pb(1))]),
        ([(pe(3), pb(2))], [(pb(3), pe(2))]),
        ([(pe(3), pb(3))], [(pe(2), pb(2)), (pe(1), pb(1))]),
        ([(pb(2), pb(3))], [(pe(2), pe(3))]),
    ]
    println("img 规范标定（VEV vs 真实 integrate）：")
    best = (Inf, -1)
    for mode in 0:5
        err = 0.0
        for (kp, ip) in probes
            M = Matrix{ComplexF64}(I, DIM, DIM)
            for (p1, p2) in kp
                M = M * (img_mono(sort([p1, p2]), mode) + Matrix{ComplexF64}(I, DIM, DIM))
            end
            Mi = Matrix{ComplexF64}(I, DIM, DIM)
            for (p1, p2) in ip
                Mi = Mi * (img_mono(sort([p1, p2]), mode) + Matrix{ComplexF64}(I, DIM, DIM))
            end
            zimg = vevc(M * Mi)
            both = vcat(sort(collect(Tuple{Int,Int}, kp)), sort(collect(Tuple{Int,Int}, ip)))
            ztrue = GTEMPO.integrate(LAT, inj(both), vac)  # 全部因子注入单链（序 = K 后 I）
            mode == 0 && @printf("  探针 K=%s I=%s : img=%g 真值=%g %s\n",
                kp, ip, zimg, real(ztrue),
                abs(zimg - real(ztrue)) > 1e-9 ? " ✗" : "✓")
            err = max(err, abs(zimg - real(ztrue)))
        end
        @printf("  mode %d : max err = %.3e\n", mode, err)
        err < best[1] && (best = (err, mode))
    end
    @printf("选中 mode = %d (err %.3e)\n", best[2], best[1])
    mode = best[2]

    # ---------- ② 同态性检验 ----------
    println("\n同态性 img(A·B) == img(A)·img(B)：")
    rng = [(Set([0, 1]), Set([2, 3])), (Set([1, 2]), Set([0, 3])), (Set([0, 5]), Set([1, 2]))]
    for (A, B) in rng
        # 多项式积（θ²=0 + reorder 符号）：用 mul_monomials 表
        monA = [(sort(collect(A)), 1.0 + 0.0im)]
        monB = [(sort(collect(B)), 1.0 + 0.0im)]
        tA = zeros(ComplexF64, 1 << L); tA[1] = 1.0
        tA = mul_monomials!(zeros(ComplexF64, 1 << L), tA, L, Monomials(monA))
        tB = zeros(ComplexF64, 1 << L); tB[1] = 1.0
        tB = mul_monomials!(zeros(ComplexF64, 1 << L), tB, L, Monomials(monB))
        tAB = mul_monomials!(zeros(ComplexF64, 1 << L), tA, L, Monomials([(sort(collect(B)), 1.0 + 0.0im)]))
        imgP(t) = begin
            M = Matrix{ComplexF64}(I, DIM, DIM)
            for idx in 1:(1 << L)
                c = t[idx]
                c == 0 && continue
                bits = [q for q in 0:(L-1) if (idx - 1) >> q & 1 == 1]
                M .+= c .* img_mono(bits, mode)
            end
            M
        end
        err = norm(imgP(tAB) - imgP(tA) * imgP(tB))
        @printf("  A=%s B=%s : |img(AB) − img(A)img(B)| = %.3e\n", collect(A), collect(B), err)
    end
end
main()
