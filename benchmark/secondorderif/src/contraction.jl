# ============================================================================
# contraction.jl — 收缩（完全复用 GTEMPO integrate）
#
# Z = ∫ dν K · IF 由 GTEMPO 的 integrate（ExactIntegrate，zipup）完成。
# 与 U1GTEMPO2 的差异：U1 的 integrate 实现 Gaussian 测度并在
# EnvironmentCache 中做双向环境缓存；基准规模下 GTEMPO 的 ExactIntegrate
# 足够，且 K/Weyl 侧无需测度适配（见 system.jl 注释）。
# ============================================================================

"""
    partition_function(lattice, mpsK, mpsIs...; alg=GTEMPO.ExactIntegrate())

归一化分母 Z = ∫ dν K·IF。
"""
partition_function(lattice, mpsK, mpsIs...; alg = GTEMPO.ExactIntegrate()) =
    integrate(lattice, mpsK, mpsIs...; alg = alg)

"""
    expectation(lattice, mpsK_num, mpsIs...; Z, alg) -> ComplexF64

归一化期望值 = 分子积分 / Z；`Z` 缺省时现场计算。
"""
function expectation(lattice, mpsK_num, mpsIs...;
                     Z::Union{Nothing,Number} = nothing, alg = GTEMPO.ExactIntegrate())
    num = integrate(lattice, mpsK_num, mpsIs...; alg = alg)
    den = isnothing(Z) ? integrate(lattice, mpsK_num, mpsIs...; alg = alg) : Z
    return num / den
end
