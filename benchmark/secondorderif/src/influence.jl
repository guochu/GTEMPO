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

"""
    influence_functional(c::WeylMixedTime, kernel; trunc=...) -> Vector{GrassmannMPS}

混合轮廓专用 IF 装配：行循环只遍历 Weyl 格顶点（τ 格 1..Nτ、rt 格 1..Nt），跳过端点顶点。
列序 = `GTEMPO._columns`（τ 列后接 ± 列）。
"""
function influence_functional(c::WeylMixedTime, kernel::HybridizationKernel;
                              trunc::TruncationScheme = truncdimcutoff(D = 128, ϵ = 1.0e-12))
    lattice = weyl_lattice(c)
    Is = Any[]
    for b in 1:c.bands
        # 混合专用列序（τ 列升序后接 ± 列升序；排除端点顶点 Nt+1/1）；
        # 列的 band 必须 = 当前带（KernelCorrView 的带 gate）。
        cols = vcat([GrassmannContourIndex(j, :τ, b, false) for j in 1:lattice.Nτ],
                    [GrassmannContourIndex(j, br, b, false)
                     for j in 1:lattice.Nt for br in (:+, :-)])
        view = KernelCorrView(kernel, b)
        gmps = vacuumstate(ComplexF64, lattice)
        for branch in (:τ, :+, :-)
            rows = branch === :τ ? (1:lattice.Nτ) : (1:lattice.Nt)
            for j in rows
                tmp = partialif_hybrid(lattice, GrassmannContourIndex(j, branch, b, true), cols, view)
                gmps = mult!(gmps, tmp, trunc = trunc)[1]
            end
        end
        push!(Is, gmps)
    end
    return Is
end

