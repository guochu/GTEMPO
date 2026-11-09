# =============================================================================
# imag_commons.jl —— 虚轴 Trotter error scaling 共享 helpers
#
# 供 toulouse/imag_scaling.jl 与 singlemodebath/imag_scaling.jl include。
# 提供：模型哈密顿量构建、连续极限 ED 参考、线性调度离散链锚（端点半门保留
# ——环合并会破坏 order-4 的 Suzuki 条件，见 imag_theory_notes.md）、CSV 输出。
#
# 依赖：NativeExact（经 ../../src/NativeExact.jl include，须先 include）。
# =============================================================================
module ImagCommons

using LinearAlgebra
using Printf

export ModelParams, build_model, continuous_ref, discrete_anchor,
       scaling_table, print_header

"""虚轴模型参数。"""
Base.@kwdef struct ModelParams
    tag::String
    WB::Float64          # bath 模能量
    AB::Float64          # 杂化强度 |V|² = AB
    BETA::Float64        # 逆温度（β = 1：T = 1）
    EPSD::Float64        # 杂质能级（tilt 参数）
end

"""单 bath 模 + 杂质（2 模，DIM=4）的算符与 Hq（不含 tilt）。"""
struct ModelOps
    dim::Int
    ops::Vector{Matrix{ComplexF64}}
    nimp::Matrix{ComplexF64}
    Hq::Matrix{ComplexF64}
end

function build_model(m::ModelParams)
    dim = 4                              # 杂质 2 × bath 2
    ops = [annih_op(dim, k) for k in 0:1]
    # 单粒子耦合矩阵（NativeExact._coupled_h 约定：h[imp,bath]=conj(V), h[bath,imp]=V；
    # 杂质能级 0 —— Toulouse εd 仅经 tilt εd·n 进入）
    V = sqrt(m.AB)
    h = [0.0+0.0im conj(V); V m.WB+0.0im]
    Hq = zeros(ComplexF64, dim, dim)
    for i in 1:2, j in 1:2
        iszero(h[i, j]) || (Hq .+= h[i, j] .* (ops[i]' * ops[j]))
    end
    return ModelOps(dim, ops, ops[1]' * ops[1], Hq)
end

function annih_op(dim::Int, mode::Int)
    op = zeros(ComplexF64, dim, dim)
    for n in 0:(dim - 1)
        (n >> mode) & 1 == 1 || continue
        sgn = isodd(count_ones(n & ((1 << mode) - 1))) ? -1 : 1
        op[(n ⊻ (1 << mode)) + 1, n + 1] = sgn
    end
    return op
end

"""连续极限参考：tr(e^{−β(Hq + εd·n)}) 与 n 期望。"""
function continuous_ref(ops::ModelOps, m::ModelParams)
    E = exp(-m.BETA .* (ops.Hq + m.EPSD .* ops.nimp))
    z = tr(E)
    return z, real(tr(ops.nimp * E) / z)
end

"""线性调度离散链锚（端点半门保留；迹循环性 = 合法环公式，保真 Suzuki 阶）：
e^{−g1 B} e^{−c1 A} e^{−g2 B} e^{−c2 A} … e^{−g_{L+1} B}，B = εd·n。"""
function discrete_anchor(ops::ModelOps, m::ModelParams, cells, links)
    Z = Matrix{ComplexF64}(I, ops.dim, ops.dim)
    Zn = Matrix{ComplexF64}(I, ops.dim, ops.dim)
    B = m.EPSD .* ops.nimp
    for i in 1:length(cells)
        Z = exp(-links[i] .* B) * Z
        Z = exp(-cells[i] .* ops.Hq) * Z
    end
    Z = exp(-links[end] .* B) * Z
    Zn = ops.nimp * Z
    return tr(Z), real(tr(Zn) / tr(Z))
end

"""阶数扫描表 → CSV。返回 (连续参考, 每 (order,N) 的 n 误差字典)。"""
function scaling_table(ops::ModelOps, m::ModelParams, schedfn, outpath,
                       orders::Vector{Int}, Ns::Vector{Int})
    zc, nc = continuous_ref(ops, m)
    open(outpath, "w") do io
        println(io, "order,N,zn_disc,n_disc,zerr_rel,nerr")
        for order in orders, N in Ns
            sched = schedfn(order, N, m.BETA / N)
            zd, nd = discrete_anchor(ops, m, sched.cells, sched.links)
            println(io, order, ",", N, ",", zd, ",", nd, ",", abs(zd - zc) / abs(zc),
                    ",", abs(nd - nc))
        end
    end
    return (zc, nc)
end

function print_header(m::ModelParams, zc, nc)
    @printf("模型 %s（WB=%.1f AB=%.2f β=%.1f εd=%.2f）\n连续极限：Z* = %.12f  n* = %.12f\n",
            m.tag, m.WB, m.AB, m.BETA, m.EPSD, zc, nc)
end

end # module
