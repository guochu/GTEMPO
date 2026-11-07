# =============================================================================
# imag_trace_probe.jl —— 判决探针：实测 integrate() 对 imag lattice 的配对测度
#
# 背景：imag_skeleton_extract.jl 已直读 Δ* = D_full（pair 匹配），但端到端
# Z/n 未过 ⇒ integrate 的 K/I 配对语义（trace pairing 与隐式 resolution）
# 与稠密 trace 测度不同。本脚本用 GTerm 注入 basis 单项式（apply! 到 vacuum，
# 全部机器精度可控），对真实 integrate() 实测 trace(mono_K ⊗ mono_I)。
#
# 生成器位（N=2，L=6，0-based）：
#   0=slice3 eta  1=slice3 bar  2=slice2 eta  3=slice2 bar  4=slice1 eta  5=slice1 bar
# 切片对（相邻融合迹候选）：(0,1) (2,3) (4,5)。
# 注入方式：mps = apply!(exp(GTerm(p1, p2, c)), vacuum)；exp(c·θp1θp2) = 1 + c·θp1θp2。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_trace_probe.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using Z2Tensors: NoTruncation

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))

const BETA = 1.0

"""vacuum MPS + 注入指定 pair 集（每对一个 GTerm exp(c·θp1θp2)）。"""
function injected_kmps(lat, pairs::Vector{Tuple{Int,Int}}; cs::Real=1.0)
    mps = GTEMPO.vacuumstate(ComplexF64, lat)
    for (p1, p2) in pairs
        mps = GTEMPO.apply!(exp(GTEMPO.GTerm(p1 + 1, p2 + 1; coeff=ComplexF64(cs))), mps)
    end
    return mps
end

function main()
    nl = 2
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    L = length(lat)
    @printf("L=%d  位图：0=s3e 1=s3b 2=s2e 3=s2b 4=s1e 5=s1b\n", L)

    Z0 = GTEMPO.integrate(lat, injected_kmps(lat, Tuple{Int,Int}[]),
                          injected_kmps(lat, Tuple{Int,Int}[]))
    @printf("探针 0  K=∅ I=∅                       : Z = %.6f（vacuum×vacuum）\n", Z0)

    probes = [
        (tag="1  K={(0,1)} I=∅",   kp=[(0, 1)], ip=Tuple{Int,Int}[], kc=1.0),
        (tag="2  K=∅ I={(0,1)}",   kp=Tuple{Int,Int}[], ip=[(0, 1)], kc=1.0),
        (tag="3  K={(0,1)} I={(0,1)}", kp=[(0, 1)], ip=[(0, 1)], kc=1.0),
        (tag="4  K={(2,3)} I={(0,1)}", kp=[(2, 3)], ip=[(0, 1)], kc=1.0),
        (tag="5  K={(0,3)} I=∅",   kp=[(0, 3)], ip=Tuple{Int,Int}[], kc=1.0),
        (tag="6  K={(0,1),(2,3)} I=∅", kp=[(0, 1), (2, 3)], ip=Tuple{Int,Int}[], kc=1.0),
        (tag="7  K={(0,1),(2,3)} I={(4,5)}", kp=[(0, 1), (2, 3)], ip=[(4, 5)], kc=1.0),
        (tag="8  K={(0,1)} I={(0,1),(2,3)}", kp=[(0, 1)], ip=[(0, 1), (2, 3)], kc=1.0),
        (tag="9  K={(0,5)} I=∅",   kp=[(0, 5)], ip=Tuple{Int,Int}[], kc=1.0),
        (tag="10 K={(1,4)} I=∅",   kp=[(1, 4)], ip=Tuple{Int,Int}[], kc=1.0),
        (tag="11 K={(0,1)}(c=-1) I=∅", kp=[(0, 1)], ip=Tuple{Int,Int}[], kc=-1.0),
        (tag="12 K={(0,1)} I={(2,3)}", kp=[(0, 1)], ip=[(2, 3)], kc=1.0),
        (tag="13 K={(4,5)} I={(4,5)}", kp=[(4, 5)], ip=[(4, 5)], kc=1.0),
        (tag="14 K={(0,2)} I={(1,3)} ⇒ r₁²", kp=[(0, 2)], ip=[(1, 3)], kc=1.0),
        (tag="15 K={(1,3)} I={(0,2)} ⇒ r₂²", kp=[(1, 3)], ip=[(0, 2)], kc=1.0),
        (tag="16 K={(0,2),(4,5)} I={(1,3)} ⇒ r₁²·2", kp=[(0, 2), (4, 5)], ip=[(1, 3)], kc=1.0),
    ]
    for pr in probes
        mk = injected_kmps(lat, pr.kp, cs=pr.kc)
        mi = injected_kmps(lat, pr.ip)
        z = GTEMPO.integrate(lat, mk, mi)
        @printf("探针 %-40s : Z = %12.6f\n", pr.tag, z)
    end
end
main()
