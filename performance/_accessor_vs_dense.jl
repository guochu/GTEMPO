# 访问器 ↔ Dense 装配逐元对照（秒级）：三轮廓 normal + BCS
# 访问器读数 = [a_{i} 行, ā_{j} 列] 规则（与 partialif 同）；
# Dense.Δ 的 2×节点布局 (a, ā) → Δ[2i, 2j-1]
using Random
using GTEMPO
using Z2Tensors
using ImpurityModelBase
using QuAPI
using QuAPI.Dense

Random.seed!(12354)

function compare(name, corr; Nτ=nothing, Nt=nothing)
	H = GrassmannHybridizationMatrix(corr)
	Δ = Matrix(H.Δ)
	legs = basis(corr)
	n = length(legs)
	maxd = 0.0
	for (p, li) in enumerate(legs), (q, lj) in enumerate(legs)
		# IF 系数四种组合（对应 _delta_terms）：(ā↑,a↑)(ā↓,a↓) 正常 +
		# (ā↑,ā↓)(a↓,a↑) 反常；(ā↑,a↓)/(ā↓,a↑) 不是 IF 项（dense=0）
		ok = (li.conj && !lj.conj && li.band == lj.band) ||                            # (ā↑,a↑)/(ā↓,a↓)
			 (li.conj && lj.conj && isodd(li.band) && !isodd(lj.band)) ||              # (ā↑,ā↓)
			 (!li.conj && !lj.conj && !isodd(li.band) && isodd(lj.band))               # (a↓,a↑)
		ok || continue
		d = abs(corr[Base.conj(li), Base.conj(lj)] - Δ[p, q])
		if d > 1e-10
			println("   [" , p, ",", q, "] (", li.branch, li.j, ",b", li.band, li.conj ? ",ā" : ",a", ")",
					"(", lj.branch, lj.j, ",b", lj.band, lj.conj ? ",ā" : ",a", ") acc=",
					round(corr[Base.conj(li), Base.conj(lj)], digits=6),
					"  dense=", round(Δ[p, q], digits=6))
		end
		maxd = max(maxd, d)
	end
	println(name, ": n=", n, "  max|访问器-dense|=", maxd)
	return maxd
end

bath = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=0.6, μ=0.7)
bcs = bcsbath(bath, Δ=0.6)

println("== normal")
compare("normal τ   (N=3)", fermionic_Δτ(semicircular(), β=1.0, N=3, μ=0.7))
compare("normal real(N=3)", fermionic_Δt(semicircular(), β=1.0, t=0.3, N=3, μ=0.7))
compare("normal mixed   ", fermionic_Δm(semicircular(), β=1.0, t=0.3, Nτ=2, Nt=1, μ=0.7))
println("== BCS")
compare("BCS τ   (N=1)", bcs_Δτ(DiscreteSpectrum([1.0],[0.25]); β=2.0, N=1, μ=0.0, Δ=0.5))
compare("BCS τ   (N=2)", bcs_Δτ(DiscreteSpectrum([1.0],[0.25]); β=2.0, N=2, μ=0.0, Δ=0.5))
compare("BCS real(N=1)", bcs_Δt(DiscreteSpectrum([1.0],[0.25]); β=2.0, t=0.1, N=1, μ=0.0, Δ=0.5))
compare("BCS mixed   ", bcs_Δm(DiscreteSpectrum([1.0],[0.25]); β=2.0, t=0.1, Nτ=1, Nt=1, μ=0.0, Δ=0.5))
