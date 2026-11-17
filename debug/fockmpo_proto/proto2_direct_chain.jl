# fockmpo 原型 v3：统一 [band1-first] 直构闭链 vs 现行 sysdynamics todense 对拍
#
# 键荷记账（全部融合恒等式已验证）：
#   外键 (流 s, 线程 t₀)：χ = q(s) − q(t₀)（四差荷——线程荷在差中抵消，全 t₀ 可存）
#   格内键 (p_j, a=s_{k−1}, b=s_k, t₀)：χ = q(b)+q(a)−q(p_j)−q(t₀)
#   真空边界 + 迹闭合：cell-1 输入流 ≡ t₀（经 ρ₀ 列指标）、cell-L 输出流 ≡ t₀（复合行指标）
#   ρ₀ 折 cell-1 的 L_B·S_in·ρ₀ 复合（a-槽哨兵 0）；cell-L 的 b-槽哨兵 0
using GTEMPO, TensorKit, LinearAlgebra
import GrassmannMPS as GM
import TensorKit as TK

const BANDS = parse(Int, (length(ARGS) >= 1 ? ARGS[1] : "1"))
const SYMNAME = length(ARGS) >= 2 ? ARGS[2] : "u1"
const N = parse(Int, length(ARGS) >= 3 ? ARGS[3] : "2")
const ORDER = parse(Int, length(ARGS) >= 4 ? ARGS[4] : "1")
const SYM = SYMNAME == "z2" ? Z2Symmetry() : SYMNAME == "u1" ? U1Symmetry() :
            SYMNAME == "u1u1" ? U1U1Symmetry() : Z2Z2Symmetry()

lattice = GrassmannLattice(N=N, δτ=0.25, bands=BANDS, contour=:imag, symmetry=SYM)
model = BANDS == 1 ? ImpurityHamiltonian([tunneling(1, 1, coeff=0.5)]; bands=1) :
                    KanamoriIM(U=1.0, J=0.2, norb=1)

K_base = sysdynamics(lattice, model; order=ORDER)
t_base = GM.todense(K_base)

G = TK.sectortype(lattice)
dim = 2^BANDS
fqs = fock_charges(G, BANDS)
z = one(G)
cells = GTEMPO.chronological_cells(lattice)
L = length(cells)

_d = GTEMPO._fock_annihilators(BANDS)
q = ORDER == 2 ? 0.5 : 1.0
locals = [[Matrix{ComplexF64}(q * I - _d[b]' * _d[b]), -_d[b]', _d[b],
           ORDER == 2 ? Matrix{ComplexF64}(q * I + _d[b]' * _d[b]) : Matrix{ComplexF64}(_d[b]' * _d[b])]
          for b in 1:BANDS]
ρ₀ = Matrix{ComplexF64}(I, dim, dim)   # imag：热迹 I

function S_factors(br)
    dt = 0.25
    S = Matrix(GTEMPO._model_propagator(model, br, ORDER == 2 ? dt / 2 : dt, nothing))
    return ORDER == 2 ? (S, S) :
           (br === :- ? (Matrix{ComplexF64}(I, dim, dim), S) : (S, Matrix{ComplexF64}(I, dim, dim)))
end

band_el(b, s, x, y) = begin
    # 全基矩阵元 ⟨x|L_b[s]|y⟩（JW 弦符号在完整 Fock 基里；带局域性 = 非 b 位必同）
    for bb in 1:BANDS
        bb == b && continue
        ((x >> (BANDS - bb)) & 1) == ((y >> (BANDS - bb)) & 1) || return 0.0 + 0.0im
    end
    locals[b][s+1][x+1, y+1]
end
_withband(a, b, v) = v == 0 ? (a & ~(1 << (BANDS - b))) : (a | (1 << (BANDS - b)))

# ---------------------------------------------------------------------------
# 键空间（四差荷 + 局部枚举）
# ---------------------------------------------------------------------------
states = 0:dim-1
qof(st) = fqs[st + 1]
qd(a, b) = GTEMPO.qsub(qof(a), qof(b))                       # q(a)−q(b)

outgroups = Dict{typeof(z),Vector{Tuple{Int,Int}}}()
for s in states, t in states
    push!(get!(outgroups, qd(s, t), Tuple{Int,Int}[]), (s, t))
end
outersect = sort!(collect(keys(outgroups)); by = GTEMPO.sortkey)
outerspace = Vect[G](d => length(outgroups[d]) for d in outersect)
outloc = Dict{Tuple{Int,Int},Tuple{typeof(z),Int}}()
for d in outersect, (i, p) in enumerate(outgroups[d])
    outloc[p] = (d, i)
end

quadgroups = Dict{typeof(z),Vector{NTuple{4,Int}}}()
for p in states, a in states, b in states, t in states
    χ = GTEMPO.qsub(GTEMPO.qsub(GTEMPO.qadd(qof(b), qof(a)), qof(p)), qof(t))
    push!(get!(quadgroups, χ, NTuple{4,Int}[]), (p, a, b, t))
end
quadsect = sort!(collect(keys(quadgroups)); by = GTEMPO.sortkey)
quadspace = Vect[G](d => length(quadgroups[d]) for d in quadsect)
quadloc = Dict{NTuple{4,Int},Tuple{typeof(z),Int}}()
for d in quadsect, (i, p) in enumerate(quadgroups[d])
    quadloc[p] = (d, i)
end

slot_tag(E, s) = begin
    b = GM.grassmann_basis(E)
    s == 0 ? b.vac : s == 3 ? b.pair : s == 1 ? b.eta : b.etabar
end
_physE(b) = GTEMPO.station_space(G, b)
koszul(slots, b) = prod((-1)^(count_ones(slots[bp]) * count_ones(slots[b])) for bp in 1:(b-1); init=1)

# ---------------------------------------------------------------------------
# 逐格站构造（[band1..B] 直出；B=1 单站兼两端）
# ---------------------------------------------------------------------------
data = Vector{Any}()
for (k, (br, jj)) in enumerate(cells)
    Sout, Sin = S_factors(br)
    first_cell = k == 1
    last_cell = k == L
    for j in 1:BANDS
        E = _physE(j)
        spaceL = j == 1 ? (first_cell ? oneunit(outerspace) : outerspace) : quadspace
        spaceR = j == BANDS ? (last_cell ? oneunit(outerspace) : outerspace) : quadspace
        v = TK.zeros(ComplexF64, spaceL ⊗ E ← spaceR)
        for m in 0:dim^2-1
            slots = [(m >> (2 * (bb - 1))) & 3 for bb in 1:BANDS]
            tag = slot_tag(E, slots[j])
            g = koszul(slots, j)
            if j == 1 && j == BANDS
                # 单站全复合 ⟨b|S_out·M·S_in·[ρ₀]|a⟩；cell-1 输入 ≡ t₀（ρ₀ 列）、cell-L 行 ≡ t₀
                M = [band_el(1, slots[1], x, y) for x in states, y in states]
                comp = Sout * M * Sin
                for t in states, a in states, b in states
                    first_cell && (a == t || continue)
                    last_cell && (b == t || continue)
                    val = 0.0 + 0.0im
                    if first_cell
                        for w in states
                            val += comp[b+1, w+1] * ρ₀[w+1, t+1]
                        end
                    else
                        val = comp[b+1, a+1]
                    end
                    iszero(val) && continue
                    (dL, iL) = first_cell ? (z, 1) : outloc[(a, t)]
                    (dR, iR) = last_cell ? (z, 1) : outloc[(b, t)]
                    v[(dL, tag.sector, TK.dual(dR))][iL, tag.slot, iR] = g * val
                end
            elseif j == 1
                # ⟨b|S_out·L₁(o₁)|p₁⟩；cell-1：a-槽 ≡ t₀（哨兵携带）；cell-L：行 ≡ t₀、b-槽哨兵
                for t in states
                    for a in (first_cell ? (0,) : states), b in (last_cell ? (t,) : states), p1 in states
                        row = last_cell ? t : b
                        val = 0.0 + 0.0im
                        for u in (0, 1)
                            val += Sout[row+1, _withband(p1, 1, u)+1] *
                                   band_el(1, slots[1], _withband(p1, 1, u), p1)
                        end
                        iszero(val) && continue
                        acarry = first_cell ? t : a
                        (dL, iL) = first_cell ? (z, 1) : outloc[(a, t)]
                        dR = GTEMPO.qsub(GTEMPO.qsub(GTEMPO.qadd(qof(row), qof(acarry)), qof(p1)), qof(t))
                        (dR, iR) = quadloc[(p1, acarry, row, t)]
                        v[(dL, tag.sector, TK.dual(dR))][iL, tag.slot, iR] = g * val
                    end
                end
            elseif j == BANDS
                # ⟨p_{B−1}|L_B·S_in·[ρ₀]|a⟩；cell-1 输入 ≡ t₀ 折 ρ₀、cell-L aR 真空
                for t in states
                    for a in (first_cell ? (t,) : states), b in (last_cell ? (t,) : states), pm in states
                        val = 0.0 + 0.0im
                        for u in (0, 1)
                            pb = _withband(pm, BANDS, u)
                            if first_cell
                                for w in states
                                    val += band_el(BANDS, slots[BANDS], pm, pb) *
                                           Sin[pb+1, w+1] * ρ₀[w+1, t+1]
                                end
                            else
                                val += band_el(BANDS, slots[BANDS], pm, pb) * Sin[pb+1, a+1]
                            end
                        end
                        iszero(val) && continue
                        (dL, iL) = quadloc[(pm, a, b, t)]
                        (dR, iR) = last_cell ? (z, 1) : outloc[(b, t)]
                        v[(dL, tag.sector, TK.dual(dR))][iL, tag.slot, iR] = g * val
                    end
                end
            else
                # ⟨p_{j−1}|L_j(o_j)|p_j⟩ 四元组直通
                for t in states, a in states, b in states, pm in states, pp in states
                    val = band_el(j, slots[j], pm, pp)
                    iszero(val) && continue
                    (dL, iL) = quadloc[(pm, a, b, t)]
                    (dR, iR) = quadloc[(pp, a, b, t)]
                    v[(dL, tag.sector, TK.dual(dR))][iL, tag.slot, iR] = g * val
                end
            end
        end
        push!(data, v)
    end
end

mps = GM.GMPS(typeof(data[1])[t for t in data])
chrono = Dict(cell => k for (k, cell) in enumerate(cells))
perm = Int[]
for cell in GTEMPO._storage_cells(lattice)
    k = chrono[cell]
    append!(perm, [(k - 1) * BANDS .+ (1:BANDS)...])
end
GM.reorder!(mps, perm; trunc=GM.truncrelerr(1e-14))

t_new = GM.todense(mps)
err = maximum(abs.(t_new .- t_base))
println("bands=$BANDS sym=$SYMNAME N=$N order=$ORDER  max |Δtodense| = ", err)
err < 1e-12 ? println("✅ 对拍通过") : begin
    println("❌ 不匹配")
    println("base: ", t_base)
    println("new : ", t_new)
end
