# Phase 0 golden fixtures dump（旧后端一次性 dump；换轨后旧栈不再可加载，
# 本脚本是权威锚点的唯一来源，必须在动任何 src 之前运行）。
#
# 运行：julia --project=. test/fixtures/dump_golden.jl
#
# golden 铁律（plan/gtempo_grassmannmps_backend_plan.md §Phase 0）：
#   - 键维序列跨栈必不同——只 dump 稠密系数/物理量（无截断算例，1e-10 档）；
#   - 只锚 Z2 通道（旧后端唯一扇区结构，与新实现 symmetry=:z2 扇区同构）；
#   - 每个 case 含：lattice 位置映射、corr 核块、K/IF 链稠密系数（小算例）、
#     Zvalue、occupation/gf 全端点序列（大算例只存物理量）。
#
# 稠密系数语义：`_gmps_dense_flat`（与 test/api/quapi_crosscheck.jl 的
# _gmps_to_dense 同源）——扁平后第 p 站轴变化最快（axis p = 链站 p、
# 1=空/2=占据），已含 scaling^L。新栈 checker 用 GM.todense + lattice
# 位置映射重排后对拍。
#
# 轮廓初态约定（与 test/normalbath 一致）：
#   虚时 = 热迹（β 由轮廓给出）；实时 = 有限 β 热初态 systhermalstate!;
#   混合 = τ 段建热态 + rt 段演化（无 systhermalstate）。

using GTEMPO
using GTEMPO: branches
using Z2Tensors: fusiontrees, dim, sectors, NoTruncation, truncdimcutoff
using GTEMPO: space_l, space_r
using ImpurityModelBase
using Random

Random.seed!(12354)

# ---------------------------------------------------------------- 极简 JSON
jsonnum(x::Float64) = isnan(x) ? "\"NaN\"" : isinf(x) ? (x > 0 ? "\"Inf\"" : "\"-Inf\"") : string(x)
jsonnum(x::Int) = string(x)
jsonnum(x::Bool) = x ? "true" : "false"
jstr(s) = "\"" * s * "\""
jarr(v::AbstractVector) = "[" * join(jsonnum.(v), ",") * "]"
jcplx(z::Complex) = "[" * jsonnum(real(z)) * "," * jsonnum(imag(z)) * "]"
jcplx(z::Real) = "[" * jsonnum(Float64(z)) * ",0.0]"
jcplxarr(v) = "[" * join(jcplx.(collect(v)), ",") * "]"
jpair(k, v) = Pair{String,String}(k, v)
jobj(items::Pair{String,<:Any}...) = "{" * join((p.first * ":" * string(p.second) for p in items), ",") * "}"
write_json(path, pairs::Pair{String,<:Any}...) = open(path, "w") do io
    write(io, jobj(pairs...))
    println("WROTE ", path, "  (", filesize(path), " bytes)")
end

# ------------------------------------------------- 旧栈稠密系数（quapi_crosscheck 同源）
_sector_offset(V, c) = sum(dim(V, d) for d in sectors(V) if d.n < c.n; init=0)

function _site_matrix(A, σ::Int)
    Vl, Vr = space_l(A), space_r(A)
    M = zeros(GTEMPO.scalartype(A), dim(Vl), dim(Vr))
    for (fl, fr) in fusiontrees(A)
        sl, sp = fl.uncoupled
        sp.n == σ || continue
        fl.coupled.n == fr.coupled.n || continue
        sr = fr.uncoupled[1]
        ol, o = _sector_offset(Vl, sl), _sector_offset(Vr, sr)
        dl, dr = dim(Vl, sl), dim(Vr, sr)
        M[ol+1:ol+dl, o+1:o+dr] .= reshape(collect(A[fl, fr]), dl, dr)
    end
    return M
end

"GMPS 稠密系数（含 scaling^L）→ 扁平复数组（站 1 轴变化最快）"
function _gmps_dense_flat(psi)
    L = length(psi)
    Ms = [[_site_matrix(psi[p], σ) for σ in 0:1] for p in 1:L]
    out = zeros(GTEMPO.scalartype(psi), ntuple(_ -> 2, L))
    for I in CartesianIndices(out)
        v = Ms[1][I[1]]
        for p in 2:L
            v = v * Ms[p][I[p]]
        end
        out[I] = v[1, 1]
    end
    return vec(out) .* (GTEMPO.scaling(psi)^L)
end

_scaling_json(x) = jobj(jpair("re", jsonnum(real(x))), jpair("im", jsonnum(imag(x))))

# ---------------------------------------------------------------- 算例装配
const OUTDIR = joinpath(@__DIR__, "port_golden")
mkpath(OUTDIR)

"(j,conj,branch,band) → 链站 的全量映射（checker 据此重排稠密系数）"
function indexmap_entries(lat)
    entries = String[]
    for b in branches(lat), band in 1:lat.bands
        kmax = b === :τ ? lat.kτ : lat.kt
        for j in 1:kmax, c in (false, true)
            push!(entries, jobj(jpair("j", jsonnum(Int(j))),
                                jpair("conj", jsonnum(Bool(c))),
                                jpair("branch", jstr(String(b))),
                                jpair("band", jsonnum(Int(band))),
                                jpair("site", jsonnum(Int(lat[GrassmannContourIndex(j, b, band, c)])))))
        end
    end
    return entries
end

"corr 核锚点：leg 序 = [τ 1..Nτ][+ 1..Nt][− 1..Nt]（band 1）；
 M_aa[i,j] = corr[aᵢ, āⱼ]、M_aa2[i,j] = corr[āᵢ, aⱼ]（QuAPI 双约定全 pin）"
function corr_matrices(corr, lat)
    legs = GrassmannContourIndex[]
    :τ in branches(lat) && append!(legs, [GrassmannContourIndex(j, :τ, 1, false) for j in 1:lat.Nτ])
    :+ in branches(lat) && append!(legs, [GrassmannContourIndex(j, b, 1, false) for j in 1:lat.kt for b in (:+, :-)])
    n = length(legs)
    conjlegs = [GrassmannContourIndex(x.j, x.branch, x.band, true) for x in legs]
    M1 = [corr[legs[i], conjlegs[j]] for i in 1:n, j in 1:n]   # a-ā
    M2 = [corr[conjlegs[i], legs[j]] for i in 1:n, j in 1:n]   # ā-a
    return n, jcplxarr(vec(M1)), jcplxarr(vec(M2))
end

"单一算例全量 dump；dense::Bool 控制是否存 K/IF 稠密系数表（仅小算例）"
function dump_case(name; contour::Symbol, lat, model, bath, alg, dense::Bool, thermal_beta::Union{Nothing,Float64})
    corr = correlationfunction(bath, lat)
    mpsI1 = hybriddynamics(lat, corr, alg)
    Is = [mpsI1]
    for b in 2:lat.bands
        push!(Is, swapband(Is[1], lat, 1, b))
    end
    mpsK = sysdynamics(lat, model, trunc=NoTruncation())
    for band in 1:lat.bands
        mpsK = boundarycondition!(mpsK, lat, band=band, trunc=NoTruncation())
    end
    if thermal_beta !== nothing
        mpsK = systhermalstate!(mpsK, lat, model, trunc=NoTruncation(), β=thermal_beta)
    end
    Z = integrate(lat, mpsK, Is...)

    fields = Pair{String,<:Any}[]
    push!(fields, jpair("kind", jstr(String(contour))))
    push!(fields, jpair("bands", jsonnum(Int(lat.bands))))
    push!(fields, jpair("alg", jstr(string(nameof(typeof(alg))))))
    push!(fields, jpair("L", jsonnum(Int(length(lat)))))
    push!(fields, jpair("thermal_beta", thermal_beta === nothing ? "null" : jsonnum(thermal_beta)))
    push!(fields, jpair("indexmap", "[" * join(indexmap_entries(lat), ",") * "]"))
    nlegs, m1, m2 = corr_matrices(corr, lat)
    push!(fields, jpair("corr_nlegs", jsonnum(nlegs)))
    push!(fields, jpair("corr_a_abar", m1))
    push!(fields, jpair("corr_abar_a", m2))
    push!(fields, jpair("Z", jcplx(Z)))
    if dense
        push!(fields, jpair("K_dense", jcplxarr(_gmps_dense_flat(mpsK))))
        push!(fields, jpair("K_scaling", _scaling_json(GTEMPO.scaling(mpsK))))
        for (b, I) in enumerate(Is)
            push!(fields, jpair("IF_dense_b$b", jcplxarr(_gmps_dense_flat(I))))
            push!(fields, jpair("IF_scaling_b$b", _scaling_json(GTEMPO.scaling(I))))
        end
    end
    # occupation 全端点 × 全带（rt 顶点 branch=:+，τ 顶点 branch=:τ）
    if :τ in branches(lat)
        for b in 1:lat.bands
            ns = [occupation(lat, i, mpsK, Is...; band=b, branch=:τ, Z=Z) for i in 1:lat.kτ]
            push!(fields, jpair("occupation_tau_b$b", jcplxarr(ns)))
        end
    end
    if :+ in branches(lat)
        for b in 1:lat.bands
            ns = [occupation(lat, i, mpsK, Is...; band=b, branch=:+, Z=Z) for i in 1:lat.kt]
            push!(fields, jpair("occupation_rt_b$b", jcplxarr(ns)))
        end
    end
    # gf 全端点序列（与 test/normalbath/util.jl 的 series 助手同约定）
    if contour === :imag
        gτ = [gf(lat, (GrassmannContourIndex(i, conj=false, branch=:τ, band=1),
                       GrassmannContourIndex(1, conj=true, branch=:τ, band=1)), mpsK, Is..., Z=Z)
              for i in 1:lat.kτ]
        push!(fields, jpair("gtau", jcplxarr(gτ)))
    elseif contour === :real
        gt = [-im * gf(lat, (GrassmannContourIndex(k, conj=false, branch=:+, band=1),
                             GrassmannContourIndex(1, conj=true, branch=:+, band=1)), mpsK, Is..., Z=Z)
              for k in 1:lat.kt]
        lt = [-im * gf(lat, (GrassmannContourIndex(1, conj=true, branch=:-, band=1),
                             GrassmannContourIndex(k, conj=false, branch=:+, band=1)), mpsK, Is..., Z=Z)
              for k in 1:lat.kt]
        push!(fields, jpair("greater", jcplxarr(gt)))
        push!(fields, jpair("lesser", jcplxarr(lt)))
    else
        gt = [-im * gf(lat, (GrassmannContourIndex(k, conj=false, branch=:+, band=1),
                             GrassmannContourIndex(1, conj=true, branch=:+, band=1)), mpsK, Is..., Z=Z)
              for k in 1:lat.kt]
        lt = [im * gf(lat, (GrassmannContourIndex(1, conj=true, branch=:-, band=1),
                            GrassmannContourIndex(k, conj=false, branch=:+, band=1)), mpsK, Is..., Z=Z)
              for k in 1:lat.kt]
        gτ = [gf(lat, (GrassmannContourIndex(k, conj=false, branch=:τ, band=1),
                       GrassmannContourIndex(1, conj=true, branch=:τ, band=1)), mpsK, Is..., Z=Z)
              for k in 1:lat.kτ]
        push!(fields, jpair("greater", jcplxarr(gt)))
        push!(fields, jpair("lesser", jcplxarr(lt)))
        push!(fields, jpair("gtau", jcplxarr(gτ)))
    end
    write_json(joinpath(OUTDIR, "port_golden_$name.json"), fields...)
end

# PartialIF 的 trunc 字段限定 TruncateDimCutoff；D=4096/ϵ=0 对 L≤16 链 = 事实无截断
partialif = PartialIF(trunc=truncdimcutoff(D=4096, ϵ=0.0))
# ExactTTIIF 短链：默认 Prony n=15 需要长数据，小算例配小 n
exactttiif_imag = ExactTTIIF(alg_expan=OverDeterminedProny(n=3, tol=1e-12, verbosity=0),
                             alg_mult=SVDCompression(trunc=NoTruncation()), verbosity=0)
exactttiif_real = ExactTTIIF(alg_expan=OverDeterminedProny(n=2, tol=1e-12, verbosity=0),
                             alg_mult=SVDCompression(trunc=NoTruncation()), verbosity=0)

# ---------------------------------------------------------------- case 矩阵
ω, α, μ, βrt = 1.0, 0.5, 0.7, 1.0
spec = DiracDelta(ω=ω, α=α)

# ---- 虚时：热迹（β 由轮廓给出） ----
let lat = GrassmannLattice(N=2, δτ=0.3, contour=:imag)
    dump_case("imag_n2_toulouse_partialif"; contour=:imag, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=true, thermal_beta=nothing)
end
let lat = GrassmannLattice(N=4, δτ=0.2, contour=:imag)
    dump_case("imag_n4_toulouse_exactttiif"; contour=:imag, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=exactttiif_imag, dense=true, thermal_beta=nothing)
end
let lat = GrassmannLattice(N=2, δτ=0.3, bands=2, contour=:imag)
    dump_case("imag_n2_anderson_partialif"; contour=:imag, lat,
              model=AndersonIM(U=0.7, ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=true, thermal_beta=nothing)
end
let lat = GrassmannLattice(N=8, δτ=0.1, contour=:imag)
    dump_case("imag_n8_toulouse_partialif"; contour=:imag, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=false, thermal_beta=nothing)
end
let lat = GrassmannLattice(N=6, δτ=0.1, bands=2, contour=:imag)
    dump_case("imag_n6_anderson_partialif"; contour=:imag, lat,
              model=AndersonIM(U=0.7, ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=false, thermal_beta=nothing)
end

# ---- 实时：有限 β 热初态 ----
let lat = GrassmannLattice(N=2, δt=0.2, contour=:real)
    dump_case("real_n2_toulouse_partialif"; contour=:real, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=βrt, μ=μ),
              alg=partialif, dense=true, thermal_beta=βrt)
end
let lat = GrassmannLattice(N=3, δt=0.15, contour=:real)
    dump_case("real_n3_toulouse_exactttiif"; contour=:real, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=βrt, μ=μ),
              alg=exactttiif_real, dense=true, thermal_beta=βrt)
end
let lat = GrassmannLattice(N=2, δt=0.2, bands=2, contour=:real)
    dump_case("real_n2_anderson_partialif"; contour=:real, lat,
              model=AndersonIM(U=0.7, ϵ_d=μ), bath=fermionicbath(spec, β=βrt, μ=μ),
              alg=partialif, dense=true, thermal_beta=βrt)
end
let lat = GrassmannLattice(N=6, δt=0.1, contour=:real)
    dump_case("real_n6_toulouse_partialif"; contour=:real, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=βrt, μ=μ),
              alg=partialif, dense=false, thermal_beta=βrt)
end
let lat = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real)
    dump_case("real_n4_anderson_partialif"; contour=:real, lat,
              model=AndersonIM(U=0.7, ϵ_d=μ), bath=fermionicbath(spec, β=βrt, μ=μ),
              alg=partialif, dense=false, thermal_beta=βrt)
end

# ---- 混合：τ 段热 + rt 段演化（无 systhermalstate） ----
let lat = GrassmannLattice(Nτ=2, δτ=0.3, Nt=1, δt=0.2, contour=:mixed)
    dump_case("mixed_small_toulouse_partialif"; contour=:mixed, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=true, thermal_beta=nothing)
end
let lat = GrassmannLattice(Nτ=2, δτ=0.3, Nt=1, δt=0.2, bands=2, contour=:mixed)
    dump_case("mixed_small_anderson_partialif"; contour=:mixed, lat,
              model=AndersonIM(U=0.7, ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=true, thermal_beta=nothing)
end
let lat = GrassmannLattice(Nτ=3, δτ=0.1, Nt=2, δt=0.1, contour=:mixed)
    dump_case("mixed_toulouse_partialif"; contour=:mixed, lat,
              model=ToulouseIM(ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=false, thermal_beta=nothing)
end
let lat = GrassmannLattice(Nτ=3, δτ=0.1, Nt=2, δt=0.1, bands=2, contour=:mixed)
    dump_case("mixed_anderson_partialif"; contour=:mixed, lat,
              model=AndersonIM(U=0.7, ϵ_d=μ), bath=fermionicbath(spec, β=lat.β, μ=μ),
              alg=partialif, dense=false, thermal_beta=nothing)
end

println("golden dump complete: ", length(readdir(OUTDIR)), " files in ", OUTDIR)
