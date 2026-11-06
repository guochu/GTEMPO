# ============================================================================
# solver.jl — Problem / Solution / solve（对齐 U1GTEMPO2 的 SIAMProblem/solve）
# ============================================================================

"""
    Problem(contour, model; bath=nothing, kernel=nothing, order=2)

单模 Anderson 基准问题：轮廓（`WeylRealTime`/`WeylImagTime`）、杂质模型
（GTEMPO `ConstImpurityHamiltonian`）与单模浴。`bath` 为 ImpurityModelBase
费米浴（如 `fermionicbath(DiracDelta(ω=ω, α=α), β=β)`），也可直接传
`kernel = hybridization(bath, contour)`。`order=1/2` 分别为一阶（anti-Wick +
Lie）与 Weyl 二阶；两者共用同一 IF。
"""
struct Problem
    contour::AbstractWeylContour
    model::Any
    kernel::HybridizationKernel
    order::Int
    lattice::AbstractGrassmannLattice
    rho0::Union{Nothing,Matrix{ComplexF64}}
end

function Problem(c::AbstractWeylContour, model; bath = nothing, kernel = nothing,
                 order::Integer = 2, rho0 = nothing)
    order in (1, 2) || throw(ArgumentError("order 必须为 1 或 2"))
    (isnothing(bath) ⊻ isnothing(kernel)) ||
        throw(ArgumentError("bath 与 kernel 恰好提供一个"))
    k = isnothing(kernel) ? hybridization(bath, c) : kernel
    return Problem(c, model, k, Int(order), weyl_lattice(c),
                   isnothing(rho0) ? nothing : Matrix{ComplexF64}(rho0))
end

"""
    Solution

求解结果：K / IF / Z 及格点。`prob` 字段保留原问题以便观测量重建 K。
"""
struct Solution
    prob::Problem
    mpsK::GrassmannMPS
    mpsIs::Vector{GrassmannMPS}
    Z::ComplexF64
end

"""
    solve(prob::Problem; iftrunc=truncdimcutoff(D=128, ϵ=1e-12)) -> Solution

构造 K（`system_factor`）与 IF（`influence_functional`，复用 GTEMPO
hybriddynamics），并给出归一化分母 Z = ∫K·IF。
"""
function solve(prob::Problem; iftrunc::TruncationScheme = truncdimcutoff(D = 128, ϵ = 1.0e-12),
               kwargs...)
    mpsK = system_factor(prob.contour, prob.model; order = prob.order,
                         coupling_norm = coupling_norm(prob.kernel), rho0 = prob.rho0,
                         kwargs...)
    mpsIs = influence_functional(prob.contour, prob.kernel; trunc = iftrunc)
    Z = partition_function(prob.lattice, mpsK, mpsIs...)
    return Solution(prob, mpsK, mpsIs, Z)
end
