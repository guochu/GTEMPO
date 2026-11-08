# =============================================================================
# imag_matrixid.jl —— Δ̃ 的矩阵闭式判决
#
# 理论（u(k) 代数封闭）：I 链 = exp(Σ Δ̃_ij c†ᵢcⱼ)，
#   Δ̃ = log( ρ_K⁻¹ · e^{Ω} )，
#   ρ_K = e^{γ₁E₂₂}e^{γ₂E₃₃}e^{E₂₃}e^{−E₃₂}e^{E₁₃}
#     （门对角 γ_s = 1−e^{−τε_s}；分辨率 p1: +E₂₃、p2: −E₃₂；BC: +E₁₃），
#   Ω = log(I + D̃ext)（D̃ = 环 Schur 嵌入切片 2..k，缝行列 = 0）。
# 判决：矩阵公式 vs seamfit2 唯一数值解（Frobenius）；扫描相似变换/变体。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_matrixid.jl
# =============================================================================
using Printf
using LinearAlgebra

const BETA = 1.0
const N = 2
const WB, AB = 1.0, 0.5
VB = [sqrt(AB) + 0.0im;;]
HB = [WB;;]

"""seamfit2 的唯一数值解（3×3 复矩阵）。"""
const DT_TARGET = ComplexF64[
    -0.482645+0.0im  -1.59029+0.813296im  -2.24338+1.1119im
    -0.080908+0.0im   0.821869-0.826748im -0.077566+0.315492im
    -1.39017-0.315083im 0.313211+2.91799im  -0.076896+1.22441im]

E(k, i, j) = (M = zeros(ComplexF64, k, k); M[i, j] = 1.0; M)

function expm(M::Matrix{ComplexF64})
    F = eigen(M)
    v = F.vectors
    return v * Diagonal(exp.(F.values)) * inv(v)
end
function logm(M::Matrix{ComplexF64})
    F = eigen(M)
    v = F.vectors
    return v * Diagonal(log.(ComplexF64.(F.values))) * inv(v)
end

function main()
    k = N + 1
    h = zeros(ComplexF64, 2, 2)                    # 单模 bath：杂质 + 1 模
    h[1, 2] = conj(VB[1]); h[2, 1] = VB[1]; h[2, 2] = HB[1, 1]
    mm = 2
    cells = fill(BETA / N, N)
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0
        end
    end
    Iset = [(r-1)*mm + 1 for r in 1:N]
    Bset = [(r-1)*mm + s for r in 1:N for s in 2:mm]
    Dschur = M[Iset, Iset] - M[Iset, Bset] * (M[Bset, Bset] \ M[Bset, Iset])
    Dt = -Dschur                                   # [[−1.0169, 1.0566],[1.0566, 1.0169]]

    # 嵌入切片 2..k（缝 = 行列 1）
    Dext = zeros(ComplexF64, k, k)
    Dext[2:k, 2:k] .= Dt
    γs = [1 - exp(-cells[g]) for g in 1:N]         # 门对角（切片 g+1）

    # ---- 变体扫描：Ω 公式 × ρ_K 因子序 × 相似变换
    omegas = [
        ("Ω = log(I+Dext)", logm(I + Dext)),
        ("Ω = Dext", Dext),
        ("Ω = log(I−Dext)", logm(I - Dext)),
        ("Ω = −log(I+Dext)", -logm(I + Dext)),
    ]
    # ρ_K⁻¹ 因子序（左到右 = 多项式逆序）
    rho_inv_variants = [
        ("标准序 e^{-E13}e^{E32}e^{-E23}e^{-γ2}e^{-γ3}",
         expm(-E(k, 1, 3)) * expm(E(k, 3, 2)) * expm(-E(k, 2, 3)) *
         expm(-γs[1] * E(k, 2, 2)) * expm(-γs[2] * E(k, 3, 3))),
        ("无 BC: e^{E32}e^{-E23}e^{-γ}e^{-γ}",
         expm(E(k, 3, 2)) * expm(-E(k, 2, 3)) *
         expm(-γs[1] * E(k, 2, 2)) * expm(-γs[2] * E(k, 3, 3))),
        ("仅 Ω", expm(zeros(ComplexF64, k, k))),
    ]
    # 相似变换候选（img 字符串 = 对角相位）
    sims = [
        ("I", I),
        ("diag(1,i,i)", Diagonal(ComplexF64[1, im, im])),
        ("diag(1,−1,−1)", Diagonal(ComplexF64[1, -1, -1])),
        ("diag(1,i,−i)", Diagonal(ComplexF64[1, im, -im])),
    ]
    best = (Inf, "", "")
    for (rname, Rinv) in rho_inv_variants
        # 反解：e^{Ω} = e^{R}·e^{Δ̃target} → Ω = log(...)
        R = inv(Rinv)
        Om = logm(expm(R) * expm(DT_TARGET))
        println("\nρ_K 变体: $rname")
        println("反解 Ω = ")
        display(round.(Om, digits=6))
        # 与 log(I+Dext) 对照
        err = norm(Om - logm(I + Dext)) / norm(logm(I + Dext))
        @printf("  |Ω − log(I+D̃ext)| 相对 = %.3e\n", err)
    end
    best = (Inf, "", "")
    for (oname, Om) in omegas, (rname, Rinv) in rho_inv_variants
        P = Rinv * expm(Om)
        dt = logm(P)
        for (sname, S) in sims
            dt2 = S * dt * inv(S)
            err = norm(dt2 - DT_TARGET) / norm(DT_TARGET)
            if err < best[1]
                best = (err, "$rname | $oname | S=$sname", dt2)
            end
        end
    end
    @printf("最优相对误差 = %.3e\n  组合: %s\n", best[1], best[2])
    if best[1] < 1e-8
        println("Δ̃（公式解）= ")
        display(round.(best[3], digits=6))
    else
        println("Δ̃(公式) = ")
        display(round.(best[3], digits=6))
        println("Δ̃(target) = ")
        display(round.(DT_TARGET, digits=6))
    end
end
main()
