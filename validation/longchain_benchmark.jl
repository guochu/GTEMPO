# GTEMPO（GrassmannMPS 后端）longchain 性能基线（plan R6）。
# 工况族对齐 U1GTEMPO2 validation/longchain_benchmark.jl 的 spirit：
#   real/imag × order 1/2 × 胞数 n ∈ sizes，reps 次重复；
# stage 语义（本栈架构）：
#   K        = sysdynamics（order=2 走 weyl2order，coupling_norm=Σα）
#   IF       = hybriddynamics(PartialIF) 逐带
#   integrate= integrate_scaled(K, Is...)
#   occ      = occupation(插入重建)
#   greater  = greater(插入重建)
# 输出 CSV 列：kind,order,n,rep,stage,seconds,alloc_bytes
# 用法：julia --project=. validation/longchain_benchmark.jl [outdir]
using GTEMPO
using ImpurityModelBase
using LinearAlgebra: Diagonal, tr
using Printf

const KINDS = (:real, :imag)
const ORDERS = (1, 2)
const SIZES = (8, 16, 32)          # 单元数（链站数 = n·bands）
const BANDS = 2
const REPS = 3
const μ = 0.7; const β0 = 1.0
const DT = 0.15; const DTAU = 0.25 / 2   # δτ 取小步长保持 IF 键维温和

function bench_case(kind, order, n, reps)
    q = exp(-β0 * μ)
    rows = Tuple{String,Int,Int,Int,String,Float64,UInt}[]
    if kind === :real
        lat = GrassmannLattice(N=n, δt=DT, bands=BANDS, contour=:real, symmetry=U1())
        ρ0 = Diagonal(ComplexF64[1.0, q, q, q * q]); ρ0 /= tr(ρ0)
    else
        lat = GrassmannLattice(N=n, δτ=DTAU, bands=BANDS, contour=:imag, symmetry=U1())
        ρ0 = nothing
    end
    model = AndersonIM(U=0.0, ϵ_d=μ)
    bath = discretefermionicbath([1.0], [0.5]; β=β0)
    cn = α_total = 1.0
    for rep in 1:reps
        # K
        t0 = time(); alloc0 = Base.gc_bytes()
        K = order == 2 ? sysdynamics_weyl2order(lat, model; coupling_norm=cn, initial_state=ρ0) :
            sysdynamics(lat, model; initial_state=ρ0)
        tK = time() - t0; allocK = Base.gc_bytes() - alloc0
        push!(rows, (string(kind), order, n, rep, "K", tK, allocK))
        # IF（逐带 PartialIF）
        t0 = time(); alloc0 = Base.gc_bytes()
        Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=cn), BANDS)
        tI = time() - t0; allocI = Base.gc_bytes() - alloc0
        push!(rows, (string(kind), order, n, rep, "IF", tI, allocI))
        # integrate
        t0 = time(); alloc0 = Base.gc_bytes()
        zs = integrate_scaled(K, Is...)
        tZ = time() - t0; allocZ = Base.gc_bytes() - alloc0
        push!(rows, (string(kind), order, n, rep, "integrate", tZ, allocZ))
        # 观测量（插入重建）：occupation + greater
        t0 = time(); alloc0 = Base.gc_bytes()
        ov = occupation(lat, model, nothing, Is...; j=1, branch=kind === :real ? :+ : :τ,
                        Z=value(zs), order=order, coupling_norm=cn, initial_state=ρ0)
        tO = time() - t0; allocO = Base.gc_bytes() - alloc0
        push!(rows, (string(kind), order, n, rep, "occupation", tO, allocO))
        if kind === :real
            kernel = hybridization(lat, bath; coupling_norm=cn)
            t0 = time(); alloc0 = Base.gc_bytes()
            gv = greater(lat, model, kernel, Is...; i=2, j=0, Z=value(zs), order=order,
                         coupling_norm=cn, initial_state=ρ0)
            tG = time() - t0; allocG = Base.gc_bytes() - alloc0
            push!(rows, (string(kind), order, n, rep, "greater", tG, allocG))
        end
    end
    return rows
end

const OUTDIR = !isempty(ARGS) ? ARGS[1] : joinpath(@__DIR__, "records", "longchain-20261116")
mkpath(OUTDIR)
const OUT = joinpath(OUTDIR, "gtempo_new.csv")
open(OUT, "w") do io
    println(io, "kind,order,n,rep,stage,seconds,alloc_bytes")
    for kind in KINDS, order in ORDERS, n in SIZES
        for (kind_s, o, nn, rep, stage, t, a) in bench_case(kind, order, n, REPS)
            println(io, "$kind_s,$o,$nn,$rep,$stage,$t,$a")
            flush(io)
        end
        println(stderr, "done $kind order=$order n=$n")
    end
end
println("written: ", OUT)
