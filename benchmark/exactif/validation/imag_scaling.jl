# =============================================================================
# imag_scaling.jl —— 虚轴 Trotter error scaling 数据（order 1/2/4 × N 扫描）
#
# 模型 a（WB=1, AB=0.5, β=1, εd=0.7）。误差 = |离散调度锚 − 连续极限 ED|：
#   离散锚 = tr(∏_p e^{−c_p Hq} e^{−g_p εd n})（GMPS 管线已验证逐位复现），
#   连续参考 = tr(e^{−β(Hq + εd n)})。
# 输出 imag_scaling.csv（plot_imag_scaling.py 出图）。
# =============================================================================
using Printf
using LinearAlgebra
using Random
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

const BETA = 1.0
const WB, AB = 1.0, 0.5
const EPSD = 0.7
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

const NT = 2                      # 杂质 + 1 bath 模
const DIM = 1 << NT
const OPS = [annih_op_f(DIM, k) for k in 0:(NT-1)]
const NIMP = OPS[1]' * OPS[1]

const HQ = begin
    h = hcoupled()
    Hq = zeros(ComplexF64, DIM, DIM)
    for i in 1:NT, j in 1:NT
        iszero(h[i, j]) || (Hq .+= h[i, j] .* (OPS[i]' * OPS[j]))
    end
    Hq
end

"""线性调度离散链（端点半门保留；迹循环性 = 合法环公式，保真 Suzuki 阶）：
e^{−g1 B} e^{−c1 A} e^{−g2 B} e^{−c2 A} … e^{−g_{L+1} B}。"""
function discrete_anchor(cells, links, eps_d)
    Z = Matrix{ComplexF64}(I, DIM, DIM); Zn = Matrix{ComplexF64}(I, DIM, DIM)
    for i in 1:length(cells)
        Z = exp(-links[i] .* (eps_d .* NIMP)) * Z
        Z = exp(-cells[i] .* HQ) * Z
    end
    Z = exp(-links[end] .* (eps_d .* NIMP)) * Z
    Zn = NIMP * Z
    return tr(Z), tr(Zn) / tr(Z)
end

"""连续极限参考：tr(e^{−β(Hq + εd n)}) 与 n 期望。"""
function continuous_ref(eps_d)
    Hfull = HQ + eps_d .* NIMP
    E = exp(-BETA .* Hfull)
    return tr(E), tr(NIMP * E) / tr(E)
end

function main()
    outdir = joinpath(GT_ROOT, "benchmark", "exactif", "validation")
    zc, nc = continuous_ref(EPSD)
    @printf("连续极限：Z* = %.12f  n* = %.12f\n", zc, nc)

    Ns = [1, 2, 4, 8, 16, 32, 64]
    file = joinpath(outdir, "imag_scaling.csv")
    open(file, "w") do io
        println(io, "order,N,zn_disc,n_disc,zerr_rel,nerr")
        for order in (1, 2, 4)
            for N in Ns
                sched = NativeExact.native_schedule(order, N, BETA / N)
                zd, nd = discrete_anchor(sched.cells, sched.links, EPSD)
                zerr = abs(zd - zc) / abs(zc)
                nerr = abs(nd - nc)
                @printf("order %d  N=%-3d : n = %12.9f（err %.3e）  zrel = %.3e\n",
                        order, N, nd, nerr, zerr)
                println(io, order, ",", N, ",", zd, ",", nd, ",", zerr, ",", nerr)
                flush(stdout)
            end
        end
    end
    println("已存 $file")
end
main()
