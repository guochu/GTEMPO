# =============================================================================
# imag_psydiag.jl —— ψ_I = ρ_I|0⟩ 的宇称/内容诊断（VEV 框架判决）
#
# τ2 = ⟨0^⊗k|img(S_K)·img(I)|0^⊗k⟩ ⇒ I 链只经 ψ_I = img(I)|0⟩ 进入。
# img(偶单项式)|0⟩ 张成偶宇称子空间（dim 2^{k-1}）。本脚本对给定
# (N, order, 调度) 解 Z-行得 ψ_I^needed，报告其偶/奇分解与残差。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_psydiag.jl <NL> <order> [cells] [rlinks]
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
    Z = Matrix{ComplexF64}(I, dim, dim)
    for p in length(cells):-1:1
        W = exp(-cells[p] .* Hq) * exp(-gate_durs[p] .* Him)
        Z = W * Z
    end
    return tr(Z)
end

function main()
    NL = parse(Int, ARGS[1])
    order = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
    if length(ARGS) >= 4
        cells = [parse(Float64, x) for x in split(ARGS[3], ",")]
        rlinks = [parse(Float64, x) for x in split(ARGS[4], ",")]
    else
        sched = NativeExact.native_schedule(order, NL, BETA / NL)
        cells = sched.cells
        rlinks = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    end
    nl = length(cells)
    LAT = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(LAT)
    K = LAT.k
    DIM = 1 << K
    EVENS = Vector{Vector{Int}}([Int[]])
    for k in 2:2:L, comb in combinations(0:(L-1), k)
        push!(EVENS, collect(comb))
    end
    @printf("NL=%d order=%d nl=%d L=%d K=%d DIM=%d evensets=%d\n", NL, order, nl, L, K, DIM, length(EVENS))

    # 切片算符
    c2 = ComplexF64[0 1; 0 0]; cd2 = ComplexF64[0 0; 1 0]
    function embed2(m2::Matrix{ComplexF64}, s::Int)
        factors = [s == j ? m2 : Matrix{ComplexF64}(I, 2, 2) for j in 1:K]
        M = factors[1]
        for j in 2:K
            M = kron(M, factors[j])
        end
        return M
    end
    CS = [embed2(c2, s) for s in 1:K]
    CDS = [embed2(cd2, s) for s in 1:K]
    VAC = let v = zeros(ComplexF64, DIM); v[1] = 1.0; v end
    slice_of(pos) = K - (pos ÷ 2)
    isbar_of(pos) = (pos % 2) == 1
    function img_mono(bits::Vector{Int})
        M = Matrix{ComplexF64}(I, DIM, DIM)
        for p in sort(bits)
            M = M * (isbar_of(p) ? CDS[slice_of(p)] : CS[slice_of(p)])
        end
        return M
    end
    vevc(M) = ComplexF64(dot(VAC, M * VAC))

    """K 链多项式（门 e^{−τε}−1 @slice g+1 + 分辨率 −1 + BC）→ img 矩阵。"""
    function rhoK_of(eps_d)
        mps = GTEMPO.vacuumstate(ComplexF64, LAT)
        for g in 1:nl
            iszero(eps_d) || begin
                sg = g + 1
                gb = LAT[GrassmannContourIndex(sg, :τ, 1, true)]
                ge = LAT[GrassmannContourIndex(sg, :τ, 1, false)]
                mps = GTEMPO.apply!(exp(GTEMPO.GTerm(gb, ge;
                                     coeff=ComplexF64(exp(-rlinks[g] * eps_d) - 1.0))), mps)
            end
            s1 = g + 1
            s2 = g == nl ? 2 : g + 2
            rb = LAT[GrassmannContourIndex(s1, :τ, 1, true)]
            re = LAT[GrassmannContourIndex(s2, :τ, 1, false)]
            mps = GTEMPO.apply!(exp(GTEMPO.GTerm(rb, re; coeff=ComplexF64(-1.0))), mps)
        end
        p1 = LAT[GrassmannContourIndex(1, :τ, 1, true)]
        p2 = LAT[GrassmannContourIndex(LAT.k, :τ, 1, false)]
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1, p2; coeff=ComplexF64(-1.0))), mps)
        # 提取多项式（单 phys 腿版）
        A = ones(ComplexF64, 1, 1)
        for j in 1:length(mps)
            t = mps[j]
            spL, domR = codomain(t), domain(t)
            dL, dR = spL[1], domR[1]
            dimL = [blockdim(dL, Z2Irrep(q)) for q in 0:1]
            dimR = [blockdim(dR, Z2Irrep(q)) for q in 0:1]
            newA = zeros(ComplexF64, dimR[1] + dimR[2], 2 * size(A, 2))
            for (f1, f2) in fusiontrees(t)
                qL = f1.uncoupled[1].n; bp = f1.uncoupled[2].n; qR = f2.uncoupled[1].n
                m = t[f1, f2]
                for α in 1:dimL[qL+1], β in 1:dimR[qR+1], cfg in 1:size(A, 2)
                    newA[(qR == 0 ? β : dimR[1] + β), cfg + (bp << (j - 1))] +=
                         m[α, 1, β] * A[(qL == 0 ? α : dimL[1] + α), cfg]
                end
            end
            A = newA
        end
        vec((GTEMPO.scaling(mps)^L) .* A)
    end

    # 校验：⟨0|img(K(0))|0⟩ 应 = 2（proto 锚点）；K(0) 多项式结构
    S0 = rhoK_of(0.0)
    rhoK0 = zeros(ComplexF64, DIM, DIM)
    for idx in 1:(1 << L)
        c = S0[idx]
        c == 0 && continue
        bits = [q for q in 0:(L-1) if (idx - 1) >> q & 1 == 1]
        rhoK0 .+= c .* img_mono(bits)
    end
    @printf("校验：⟨0|img(K(0))|0⟩ = %.6f（应=2）\n", real(vevc(rhoK0)))
    println("K(0) 多项式非零项：")
    for idx in 1:(1 << L)
        abs(S0[idx]) > 1e-9 || continue
        bits = [q for q in 0:(L-1) if (idx - 1) >> q & 1 == 1]
        println("   ", bits, " : ", round(real(S0[idx]), digits=6))
    end

    # Z 行（εd 网格）
    hgrid = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]
    rows = Vector{Vector{ComplexF64}}(); rhs = ComplexF64[]
    for ed in hgrid
        S = rhoK_of(ed)
        rhoK = zeros(ComplexF64, DIM, DIM)
        for idx in 1:(1 << L)
            c = S[idx]
            c == 0 && continue
            bits = [q for q in 0:(L-1) if (idx - 1) >> q & 1 == 1]
            rhoK .+= c .* img_mono(bits)
        end
        # 行向量：r[i] = ⟨0|ρK·E_i|0⟩ 的有效形式 = ⟨0|ρK|⟩作用于 ψ：r_ψ = ⟨0|ρK|ψ⟩
        # 方程：⟨0|ρK·ψ⟩ = z。线性泛函系数 = (ρK)_{0,:}（⟨0|ρK|ψ⟩ = Σ (ρK)_{0j} ψ_j）
        push!(rows, vec(rhoK[K == 0 ? 1 : 1, :]))           # ⟨0| = 第一行（VAC[1]=1）
        push!(rhs, fock_anchor(cells, rlinks, ed)[1])
    end
    Amat = hcat(rows...)                                    # DIM × nrows
    ψ = pinv(Amat) * rhs
    res = norm(Amat * ψ - rhs)
    @printf("Z-行解 ψ：残差 = %.3e\n", res)
    # 宇称分解：偶宇称基 = img(偶集)|0⟩ 的 span —— 直接看 ψ 的 Fock 宇称
    println("ψ 非零分量（Fock 基 idx 1-based，宇称 = popcount(idx-1) mod 2）：")
    for i in 1:DIM
        abs(ψ[i]) > 1e-10 || continue
        par = count_ones(i - 1) % 2
        @printf("   |%s⟩ par=%d : %10.6f\n", bits_str(i - 1, K), par, real(ψ[i]))
    end
    even = sum(abs(ψ[i]) for i in 1:DIM if count_ones(i - 1) % 2 == 0)
    odd = sum(abs(ψ[i]) for i in 1:DIM if count_ones(i - 1) % 2 == 1)
    @printf("偶范数 = %.4f  奇范数 = %.4f\n", even, odd)
end
function bits_str(v::Int, K::Int)
    return join([string((v >> (K - j)) & 1) for j in 1:K], "")
end
main()
