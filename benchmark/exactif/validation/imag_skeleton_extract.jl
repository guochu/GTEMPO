# =============================================================================
# imag_skeleton_extract.jl —— 虚轴 GMPS 闭环最后一块：骨架稠密提取 + Δ 直读
#
# 权威入口（记忆 2026-11-07 定稿）：标定 v2 三重排除（Z 归一无辜 / gτ 读出无辜 /
# Δ 的 Schur 代数无辜）⇒ 结唯一落在管线骨架隐式测度。本脚本不猜旋钮：
#
#   E0  提取器 mps_to_dense：GMPS site 链 → canonical 稠密系数表。
#       约定（锚点标定）：poly = scaling(mps)^L × plain 链收缩；site j 载生成器
#       pos j；bond 态 flat = (sector q, 块内 α)。锚：vacuum→{∅:1}；
#       εd=0 骨架（BC=exp(GTerm(1̄,k,−1))）→ {∅:1, {1,k̄}: +1} = 1+θ₁θ₆ ✓。
#   E1  IF 基探针：Δ=E_ij 喂 hybriddynamics → 提取 pair 系数
#       （行 = corr 节点 i ↔ 切片 i+1 bar；列 = 节点 j ↔ 切片 j+1 eta），
#       线性映射 M[(i,j)] = (pos 对, 符号)——siteflip/位点序一次性读出。
#   E2  稠密目标 T_imp：imag_dense 全 kernel（εd=0，W=C）→ integrate_out(bath)
#       → 杂质扇区 canonical 表。ν_r = T_imp[∅]（bath 积分标量）。
#   E3  匹配直读：pair 系数线性（Gaussian 的 pair 系数无交叉污染——两个 bilinear
#       单项式之积 ≥ 4 腿或幂零）：
#           pairs(S) + pairs(I(Δ*)) = pairs(T_imp)/ν_r
#       ⇒ Δ*_ij = (pairs_T[pair(i,j)]/ν_r − pairs_S[pair(i,j)]) / M_sign。
#       全多项式强检验：S×I(Δ*) == ν_r·T_imp（含 4 腿系数，mul_monomials!）。
#   E4  端到端：Δ* 喂 run_pipeline（εd=0.7），Z/n vs Fock 环锚（机器精度判据）。
#       ν = 1/ν_r（εd=0 处一次锁定；S(εd)×I(Δ*) = ν_r·T_imp(εd) 应对全 εd 成立
#       ——严格 IF 的结构性主张）。
#
# 运行：julia --project=<GTEMPO 根> benchmark/exactif/validation/imag_skeleton_extract.jl
# =============================================================================
using Printf
using LinearAlgebra
using GTEMPO
using QuAPI: AbstractCorrelationFunction
using Z2Tensors: NoTruncation
using Z2Tensors: Z2Irrep, blockdim, codomain, domain, fusiontrees

const GT_ROOT = dirname(dirname(pathof(GTEMPO)))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "formb_dense.jl"))
include(joinpath(GT_ROOT, "benchmark", "exactif", "validation", "imag_dense.jl"))  # 稠密基准（include 安全）
include(joinpath(GT_ROOT, "benchmark", "exactif", "src", "NativeExact.jl"))
using .NativeExact

# imag_dense 提供：BETA HB VB EPS_D hcoupled interval_operators dense_ring fock_ring ed_n
# formb_dense 提供：Monomials mul_monomials! integrate_out kernel_monomials

# ------------------------------------------------- 环 M-Schur（正确 Iset）--
"""未拆分环 M 的 Schur 消元（杂质模 = slot s=1）。返回 (Δ = −D_full, ν_a)。"""
function ring_delta(cells::Vector{Float64})
    N = length(cells)
    h = hcoupled()
    mm = size(h, 1)
    Us = [exp(-cells[r] .* h) for r in 1:N]
    M = zeros(ComplexF64, mm * N, mm * N)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0                 # 反周期缝（M[N,N] = +U_N）
        for s2 in 1:mm, s1 in 1:mm
            M[(r-1)*mm + s2, (r-1)*mm + s1] += sgn * Us[r][s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm
            M[(r-1)*mm + s, (nxt-1)*mm + s] += 1.0  # 迹 resolution（前向对角）
        end
    end
    Iset = [(r-1)*mm + 1 for r in 1:N]              # 杂质模 = slot 1
    Bset = [(r-1)*mm + s for r in 1:N for s in 2:mm]
    D = M[Iset, Iset] - M[Iset, Bset] * (M[Bset, Bset] \ M[Bset, Iset])
    Mb = zeros(ComplexF64, (mm-1)*N, (mm-1)*N)
    hbb = exp(-(BETA / N) .* HB)
    for r in 1:N
        sgn = r == N ? +1.0 : -1.0
        for s2 in 1:mm-1, s1 in 1:mm-1
            Mb[(r-1)*(mm-1) + s2, (r-1)*(mm-1) + s1] += sgn * hbb[s2, s1]
        end
        nxt = r == N ? 1 : r + 1
        for s in 1:mm-1
            Mb[(r-1)*(mm-1) + s, (nxt-1)*(mm-1) + s] += 1.0
        end
    end
    return -D, det(M) / det(Mb), det(Mb)
end

# ------------------------------------------------------- ImagDelta 适配器 --
mutable struct ImagDelta <: AbstractCorrelationFunction
    Δ::Matrix{ComplexF64}
end
GTEMPO.TK.scalartype(::Type{<:ImagDelta}) = ComplexF64
Base.size(x::ImagDelta, args...) = size(x.Δ, args...)
function Base.getindex(x::ImagDelta, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.branch == :τ && j.branch == :τ) || throw(ArgumentError("branch must be :τ"))
    return x.Δ[i.j, j.j]                            # 符号全部进 Δ 本体（无 siteflip 旋钮）
end

# ------------------------------------------------------------ 提取器 E0 --
"""GMPS → canonical 稠密系数表（长度 2^L；位 p = 生成器 pos p+1，1-based）。

poly = scaling^L × plain 链收缩；site j 载 pos j；bond 态 flat = (q, 块内 α)。"""
function mps_to_dense(x::GTEMPO.GrassmannMPS, L::Int)
    n = length(x)
    n == L || throw(ArgumentError("site 数应等于 L=2·bands·k（每 site 一条 phys 腿），n=$n L=$L"))
    A = ones(ComplexF64, 1, 1)                      # (bond 态 × config)
    for j in 1:n
        t = x[j]
        spL, domR = codomain(t), domain(t)
        (length(spL) == 2 && length(domR) == 1) ||
            error("site 张量应为 cod=(bond,phys) dom=(bond)：site $j cod=$spL")
        dL, dR = spL[1], domR[1]
        dimL = [blockdim(dL, Z2Irrep(q)) for q in 0:1]
        dimR = [blockdim(dR, Z2Irrep(q)) for q in 0:1]
        offR = (0, dimR[1])
        newA = zeros(ComplexF64, dimR[1] + dimR[2], 2 * size(A, 2))
        for (f1, f2) in fusiontrees(t)
            qL = f1.uncoupled[1].n; bp = f1.uncoupled[2].n; qR = f2.uncoupled[1].n
            m = t[f1, f2]
            for α in 1:dimL[qL+1], β in 1:dimR[qR+1], cfg in 1:size(A, 2)
                newA[offR[qR+1] + β, cfg + (bp << (j - 1))] +=
                     m[α, 1, β] * A[α + (qL == 1 ? dimL[1] : 0), cfg]
            end
        end
        A = newA
    end
    return vec((scaling(x)^L) .* A)
end

"""稠密表 → pair 系数字典（键 = 排序 pos 对 (0-based)）。仅检查高斯性（无 >4 腿）。"""
function pair_coeffs(dense::Vector{ComplexF64}, L::Int)
    out = Dict{Tuple{Int,Int}, ComplexF64}()
    for idx in 1:length(dense)
        c = dense[idx]
        iszero(c) && continue
        bits = Int[p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
        if length(bits) == 2
            out[(bits[1], bits[2])] = get(out, (bits[1], bits[2]), 0.0) + c
        elseif length(bits) > 4
            error("发现 $(length(bits)) 腿单项式——非高斯，匹配方案失效")
        end
    end
    return out
end

"""稠密表 → Monomials（全部非零项含恒等——mul_monomials! 防全灭）。"""
function dense_to_monos(dense::Vector{ComplexF64}, L::Int)
    monos = Tuple{Vector{Int}, ComplexF64}[]
    for idx in 1:length(dense)
        c = dense[idx]
        iszero(c) && continue
        bits = Int[p for p in 0:(L-1) if (idx - 1) >> p & 1 == 1]
        push!(monos, (bits, c))
    end
    return Monomials(monos)
end

# ------------------------------------------------------------- 管线运行 --
"""GMPS τ 环管线（与 imag_pipeline_calib.jl 同构；siteflip 已并入 Δ）。

返回 (Z_pipe = integrate·ν, gτ 变体组)。"""
function run_pipeline(Δ, ν, eps_d; N=2, order=2)
    sched = NativeExact.native_schedule(order, N, BETA / N)
    links = sched.links
    ring_links = vcat(links[1] + links[end], links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(Δ); band=1, trunc=NoTruncation())
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    model = ToulouseIM(ϵ_d=eps_d)
    for g in 1:nl
        fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
        sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
        mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    Z = GTEMPO.integrate(lat, mpsK, mpsI) * ν
    b1 = GrassmannContourIndex(1, conj=true, branch=:τ, band=1)
    g_end = real(GTEMPO.gf(lat, (GrassmannContourIndex(lat.Nτ, conj=false, branch=:τ, band=1), b1),
                           mpsK, mpsI, Z=Z))
    g_1 = real(GTEMPO.gf(lat, (GrassmannContourIndex(1, conj=false, branch=:τ, band=1), b1),
                         mpsK, mpsI, Z=Z))
    return Z, (g = g_end, g1 = g_1, one_minus = 1 - g_1, one_minus_end = 1 - g_end)
end

"""构建骨架 mpsK（门 + BC，无 cc）。eps_d=0 时跳过门（恒等；且 sparse-mult
会把 phys 腿融合成 4 维，破坏 mps_to_dense 的单腿提取——内容不变）。"""
function build_skeleton(N::Int, eps_d, order::Int)
    sched = NativeExact.native_schedule(order, N, BETA / N)
    ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
    nl = length(ring_links)
    lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
    mpsK = GTEMPO.vacuumstate(ComplexF64, lat)
    if eps_d != 0
        model = ToulouseIM(ϵ_d=eps_d)
        for g in 1:nl
            fm = GTEMPO.fock_propagator(model, :τ, ring_links[g])
            sparse = GTEMPO._propagator_sparsegmps(lat, fm, g, :τ, nothing, nothing)
            mpsK = GTEMPO.mult!(mpsK, sparse, trunc=NoTruncation())[1]
        end
    end
    mpsK = GTEMPO.boundarycondition!(mpsK, lat; band=1, trunc=NoTruncation())
    return lat, mpsK
end

# -------------------------------------------------- 稠密目标表（E2）--
"""imag_dense.dense_ring 的表版本：返回积分前的全 kernel 稠密表（逐字移植）。

约定同 imag_dense（W_p kernel + 迹 resolution + 缝位反周期翻号）。"""
function dense_ring_table(Ws)
    N = length(Ws)
    nm = size(Ws[1], 1) |> x -> round(Int, log2(x))
    L = 2nm * N
    L > 22 && error("稠密规模 2^$L 过大（内存保护）")
    dense = zeros(ComplexF64, 1 << L)
    dense[1] = 1.0
    bar(p, s) = 2 * ((p - 1) * nm + s - 1)
    eta(p, s) = bar(p, s) + 1
    for p in 1:N
        kmonos = kernel_monomials([bar(p, s) for s in 1:nm],
                                  [eta(p, s) for s in 1:nm], Ws[p])
        if p == N
            kmonos = [(loc, isodd(count(b -> b in [bar(N, s2) for s2 in 1:nm], loc)) ?
                       -c : c) for (loc, c) in kmonos]
        end
        dense = mul_monomials!(zeros(ComplexF64, 1 << L), dense, L, Monomials(kmonos))
        nxt = p == N ? 1 : p + 1
        for s in 1:nm
            rmonos = Monomials([(Int[], 1.0 + 0.0im),
                                ([bar(p, s), eta(nxt, s)], -1.0 + 0.0im)])
            dense = mul_monomials!(zeros(ComplexF64, 1 << L), dense, L, rmonos)
        end
    end
    return dense, L
end

# ================================================================== 主程序 --
function main()
    NMAX = 3
    for N in 2:NMAX
        cells = fill(BETA / N, N); durs = copy(cells)
        sched = NativeExact.native_schedule(2, N, BETA / N)
        ring_links = vcat(sched.links[1] + sched.links[end], sched.links[2:end-1])
        nl = length(ring_links)
        lat = GTEMPO.GrassmannLattice(N=nl, δτ=BETA / nl, contour=:imag, bands=1)
        L = length(lat)
        println("─"^78)
        @printf("N=%d  nl=%d  L=%d\n", N, nl, L)

        # ---- E0：骨架（εd=0 门恒等 → 纯 BC 测度）
        _, mpsK0 = build_skeleton(N, 0.0, 2)
        S = mps_to_dense(mpsK0, L)
        pS = pair_coeffs(S, L)
        println("骨架 S（εd=0）：∅ = ", round(real(S[1]), digits=12),
                "  pair 项：", Dict(sort(collect(keys(pS))) .=> [round(real(pS[k]), digits=6) for k in sort(collect(keys(pS)))]))

        # ---- E1：IF 基探针（线性映射 M）
        Nτ = lat.Nτ
        Mmap = Dict{Tuple{Int,Int}, Tuple{Tuple{Int,Int}, ComplexF64}}()   # (i,j) → (pair, sign)
        for i in 1:Nτ, j in 1:Nτ
            E = zeros(ComplexF64, Nτ, Nτ); E[i, j] = 1.0
            mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(E); band=1, trunc=NoTruncation())
            pI = pair_coeffs(mps_to_dense(mpsI, L), L)
            filter!(kv -> abs(kv[2]) > 1e-10, pI)
            length(pI) == 1 || error("E[$i,$j] 探针非单 pair：$pI")
            (pr, sg) = only(pI)
            Mmap[(i, j)] = (pr, sg)
        end
        println("IF 基映射 M：")
        for i in 1:Nτ, j in 1:Nτ
            pr, sg = Mmap[(i, j)]
            @printf("  Δ[%d,%d] → pair %s 符号 %+s（切片 %d bar / 切片 %d eta）\n",
                    i, j, pr, real(sg) > 0 ? "+" : "-", i + 1, j + 1)
        end

        # ---- E2：稠密目标 T_imp（εd=0，W=C；bath 积掉）
        Ws, _ = interval_operators(cells, ring_links, 0.0)
        ztot, Ld = dense_ring_table(Ws)
        bath_gens = Int[]                                    # s=2 mode（nm=2）
        for p in 1:N
            push!(bath_gens, 4 * (p - 1) + 2, 4 * (p - 1) + 3)
        end
        Timp_full = integrate_out(ztot, Ld, bath_gens)
        ν_r = Timp_full[1]
        # 杂质位 → 管线位（ring (p,bar)=4(p−1) → 切片 p+1 bar；(p,eta)=4(p−1)+1 → 切片 p+1 eta）
        rpos2ppos(p, isbar) = lat[GrassmannContourIndex(p + 1, :τ, 1, isbar)] - 1   # 0-based
        pT = Dict{Tuple{Int,Int}, ComplexF64}()
        for (k12, c) in pair_coeffs(Timp_full, Ld)
            # 检查 pair 是否纯杂质位
            all(b -> b % 4 < 2, collect(k12)) || error("T_imp pair 含 bath 位：$k12")
            pb = rpos2ppos(k12[1] ÷ 4 + 1, k12[1] % 4 == 0)
            pe = rpos2ppos(k12[2] ÷ 4 + 1, k12[2] % 4 == 0)
            pT[tuple(sort([pb, pe])...)] = c
        end
        println("T_imp（bath 积掉）：ν_r = ", round(real(ν_r), digits=12),
                "  pair 项（管线位）：")
        for k in sort(collect(keys(pT)))
            @printf("    %s = %s\n", k, round(real(pT[k]), digits=8))
        end

        # ---- E3：直读 Δ* = (pairs_T/ν_r − pairs_S) / M_sign
        Nτ2 = Nτ
        Δstar = zeros(ComplexF64, Nτ2, Nτ2)
        for i in 1:Nτ2, j in 1:Nτ2
            pr, sg = Mmap[(i, j)]
            tval = get(pT, pr, 0.0) / ν_r
            sval = get(pS, pr, 0.0)
            Δstar[i, j] = (tval - sval) / sg
        end
        println("Δ*（直读）= ")
        display(round.(Δstar, digits=8))
        Δfull, ν_a, _ = ring_delta(cells)
        println("对照 −D_full（Schur）= ")
        display(round.(Δfull, digits=8))

        # ---- E3b：全多项式强检验 S×I(Δ*) == ν_r·T_imp（4 腿系数）
        mpsI = GTEMPO.hybriddynamics(lat, ImagDelta(Δstar); band=1, trunc=NoTruncation())
        Istar = mps_to_dense(mpsI, L)
        comp = mul_monomials!(zeros(ComplexF64, 1 << L), S, L, dense_to_monos(Istar, L))
        # ν_r·T_imp 折到管线位（高位 bath 部分恒零，直接全表对比）
        Tpipe = zeros(ComplexF64, 1 << L)
        for idx in 1:length(Timp_full)
            c = Timp_full[idx]
            iszero(c) && continue
            bits = Int[p for p in 0:(Ld-1) if (idx - 1) >> p & 1 == 1]
            all(b -> b % 4 < 2, bits) || error("T_imp 残留 bath 位：$bits")
            newbits = Int[rpos2ppos(b ÷ 4 + 1, b % 4 == 0) for b in bits]
            Tpipe[1 + sum(1 << b for b in newbits; init=0)] += c
        end
        err_max = maximum(abs.(comp .- ν_r .* Tpipe))
        nz = count(!iszero, comp) + count(!iszero, Tpipe)
        println("全多项式强检验：max |S×I(Δ*) − ν_r·T_imp| = ", @sprintf("%.3e", err_max),
                "  （非零槽位参考数 ", nz, "）")

        # ---- E4：端到端 εd=0.7（Z/n vs Fock 环锚；ν = 1/ν_r）
        eps_d = EPS_D
        FZ, FZn, _ = fock_ring(cells, ring_links, eps_d)
        zref, nref = real(tr(FZ)), real(tr(FZn) / tr(FZ))
        for νt in (1 / ν_r, ν_a, 1.0)
            z, r = run_pipeline(Δstar, νt, eps_d; N=N)
            @printf("ν=%-8s : Z = %12.8f (err %9.2e) | gτ(end)=%9.6f gτ(1)=%9.6f 1−gτ(1)=%9.6f 1−gτ(end)=%9.6f | 锚 n=%.8f\n",
                    @sprintf("%.4g", real(νt)), real(z), abs(real(z) - zref),
                    r.g, r.g1, r.one_minus, r.one_minus_end, nref)
        end
        @printf("锚：Z = %.10f  n = %.10f\n", zref, nref)
    end
end
main()
