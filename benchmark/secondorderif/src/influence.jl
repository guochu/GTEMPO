# ============================================================================
# influence.jl — 浴影响泛函（IF）
#
# IF 与一阶 GTEMPO 完全相同（物理浴、相同的 N 个全宽源单元、相同的变量约定），
# 因此整个构造复用 GTEMPO 的 hybriddynamics（PartialIF 逐行 MPO 乘法），唯一
# 的新东西是 kernels.jl 里单元级的核视图 KernelCorrView。
#
# 多带情形（如 AndersonIM bands=2，每带一个单模浴）：先按 band=1 构造，再用
# GTEMPO 的 swapband 复制到其余带（与 GTEMPO 测试 fermionic_setup 的做法一致）。
# ============================================================================

"""
    influence_functional(c::AbstractWeylContour, kernel::HybridizationKernel;
                         trunc=truncdimcutoff(D=128, ϵ=1e-12)) -> Vector{GrassmannMPS}

在 Weyl 轮廓上构造 IF = exp(Σ ā Δ a) 的 GMPS（复用 GTEMPO hybriddynamics）。
每带直接经 `KernelCorrView(kernel, b)` 独立构造（2026-11-05：swapband 复制
路径的 DP 探针给出零真空振幅，弃用；恒等浴下逐带构造与 U1 的
per-spin influence + fit_spin 语义等价且无复制步骤）。
"""
function influence_functional(c::AbstractWeylContour, kernel::HybridizationKernel;
                              trunc::TruncationScheme = truncdimcutoff(D = 128, ϵ = 1.0e-12))
    lattice = weyl_lattice(c)
    Is = [hybriddynamics(lattice, KernelCorrView(kernel, b); band = b, trunc = trunc)
          for b in 1:c.bands]
    return Is
end
