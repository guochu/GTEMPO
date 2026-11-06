# ============================================================================
# diagnostics.jl — 诊断与收敛阶拟合
# ============================================================================

"GTEMPO GMPS 的最大键维数"
maxbond(mps::GrassmannMPS) = maximum(GTEMPO.space_r.(mps.data); init = nothing) === nothing ? 1 :
                             maximum(dim(GTEMPO.space_r(t)) for t in mps.data)

"收敛阶：log-log 斜率（对齐 secondorder/convergence.jl 的拟合）"
function fit_order(hs::AbstractVector, errs::AbstractVector)
    length(hs) == length(errs) >= 2 || throw(ArgumentError("至少两个数据点"))
    x = log.(hs); y = log.(errs)
    n = length(hs)
    return (n * dot(x, y) - sum(x) * sum(y)) / (n * dot(x, x) - sum(x)^2)
end

"""
    convergence_report(prob_builder, hs, obs; reference) -> (orders, table)

`prob_builder(h)` 为给定步长构造 `Problem`，`obs(sol)` 返回观测量向量，
`reference(h)` 为同网格的精确参考。返回各步长的最大相对误差与拟合阶数。
"""
function convergence_report(prob_builder, hs, obs; reference)
    errs = Float64[]
    table = Any[]
    for h in hs
        sol = solve(prob_builder(h))
        o = obs(sol)
        ref = reference(h)
        e = norm(o - ref) / norm(ref)
        push!(errs, e)
        push!(table, (h = h, err = e, Z = sol.Z, maxbond_K = maxbond(sol.mpsK)))
    end
    return fit_order(collect(float.(hs)), errs), table
end
