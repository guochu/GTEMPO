# ============================================================================
# diag_tail_investigation.jl — 内部格结构尾部的定位判决（第十二轮）
#
# 结论：τ 两点 ~1e-3 结构尾部 = (交错 Weyl-K ⊗ IF₁ ⊗ IF₂) 分级收缩的**值**
# 与 M=1 乘积之差，绑定 IF₂ 的物理内容，一切表示变换不变。五条路线穷尽：
#   ① K 侧符号：14 特征 GF(2) 邻域穷举 + 32 结构变体（4 取向项 × 门控）
#     + 解析 χ₀ 形式（bare-joint=楔积×(−1)^χ₀ 的逐格键化）——全部无命中
#     或过矫正（χ₀ 形式残差反号）。
#   ② IF 构造：逐带直构 vs swapband 复制（上游模式）——逐位相同。
#   ③ 收缩系统序：integrate(K,IF₁,IF₂) vs integrate(K,IF₂,IF₁)——逐位相同。
#   ④ 格点置换：changeordering(BandMajorTime)（GTEMPO 精确 swap 网络）——
#     逐位相同（收缩的值对顺序不变，范畴论必然）。
#   ⑤ IF₂ 解耦：IF₂→真空 MPS ⟹ **M2−M1 = 1.5e-15 机器精确**（且与补偿
#     设置无关——IF₂ 真空使全部 c₂ 项死掉，纯 F 项恰好正确）。
# 判定：失真 = 跨带分级配对的 bilinear 结构（K-slots × IF₂-selections），
# K-only 符号在数学上不可能精确修复；上游 GTEMPO 一阶形式论同量 clean
# （sysdynamics K，relerr 2.9e-4 = 纯离散化），但其 K 是另一形式论（非
# Weyl 联合 automaton），无法移植。修复路径（未来）= GTEMPO 分级融合约定
# 的多系统 surgery，或二阶 K 的带顺序重构（Trotter 化 U 交互单元）。
# 判据工具（本脚本 + diag_j0_factor.jl）永久保留。
# 用法：julia --project=<repo> <repo>/benchmark/secondorderif/validation/diag_tail_investigation.jl <repo>
# ============================================================================
include(joinpath(ARGS[1], "benchmark", "secondorderif", "src", "SecondOrderWeyl.jl"))
using .SecondOrderWeyl
using .SecondOrderWeyl: matsubara
include(joinpath(ARGS[1], "benchmark", "secondorderif", "validation", "ed.jl"))
using .SingleModeED
using GTEMPO
using Printf

const μ, ω, α, β = 0.7, 1.0, 0.5, 1.0
const N = 4
bathA = fermionicbath(DiracDelta(ω = ω, α = α), β = β)
bath0 = fermionicbath(DiracDelta(ω = ω, α = 1e-14), β = β)
c1 = WeylImagTime(Nτ = N, dτ = β / N, bands = 1)
c2 = WeylImagTime(Nτ = N, dτ = β / N, bands = 2)
mdl2 = AndersonIM(ϵ_d = μ, U = 0.0)
sol1 = solve(Problem(c1, ToulouseIM(ϵ_d = μ); bath = bathA, order = 2))
m1 = [real(matsubara(sol1, i, 0)) for i in 1:N]

lat = weyl_lattice(c2)
dop = SecondOrderWeyl._fock_ops(2)
kernA = hybridization(bathA, c2)
kern0 = hybridization(bath0, c2)
IF_A1 = SecondOrderWeyl.influence_functional(c2, kernA)[1]
IF_A2 = SecondOrderWeyl.influence_functional(c2, kernA)[2]
IF_02 = SecondOrderWeyl.influence_functional(c2, kern0)[2]
IF2_swap = GTEMPO.swapband(IF_A1, lat, 1, 2, trunc = GTEMPO.truncdimcutoff(D = 128, ϵ = 1e-12))

K0 = system_factor(c2, mdl2; order = 2, coupling_norm = α)
Kins = [system_factor(c2, mdl2; order = 2, coupling_norm = α,
                      insertions = [WeylInsertion(:τ, i, dop[1]), WeylInsertion(:τ, 0, dop[1]')])
        for i in 1:N]

println("== U=0 判据（M2 vs M1，直接列 d@i,d†@0，N=4）==")
for (lbl, IF2) in (("IF₂=full   ", IF_A2), ("IF₂=swapband", IF2_swap), ("IF₂≈0      ", IF_02))
    Z = SecondOrderWeyl.partition_function(lat, K0, IF_A1, IF2)
    worst = 0.0
    for i in 1:N
        num = SecondOrderWeyl.partition_function(lat, Kins[i], IF_A1, IF2)
        worst = max(worst, abs(real(num / Z) - m1[i]))
    end
    @printf("  %s: worst = %.3e\n", lbl, worst)
end
println("判定：IF₂=full/swapband 残差 ~2e-3（结构尾部）；IF₂≈0 机器精确 ⟹")
println("失真绑定 IF₂ 物理内容（跨带分级配对 bilinear），K 侧不可修。")
println("diag_tail_investigation 完成")
