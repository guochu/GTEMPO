@testset "Few-mode bath, mixed time: orderings" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)

	# ED references with the interacting thermal state (built by the τ leg)
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
	gτ_ref = gτ_ed(H, a, adag, 0:δτ:β, β)

	spec = DiracDelta(ω=ω, α=α)
	bath = fermionicbath(spec, β=β)
	model = AndersonIM(U=0, μ=μ)
	for ordering in mixed_orderings
		lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, ordering=ordering)
		corr = correlationfunction(bath, lat)
		mpsI = hybriddynamics(lat, corr, trunc=trunc)
		mpsK = sysdynamics(lat, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lat)
		gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
		@test relerr(gt, gt_ed) < 3.0e-2
		@test relerr(lt, lt_ed) < 3.0e-2
		@test relerr(gτ, gτ_ref) < 1.0e-2
	end
end

@testset "Few-mode bath, mixed time: two bands with U" begin
	μ = 0.7; U = 1.0; ω = 1.0; α = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)

	# one bath mode per band
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
	gτ_ref = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, bands=2)
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, μ=μ), trunc)
	gt, lt, gτ = gtltgτ_series(lat, mpsK, Is...)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
	@test relerr(gτ, gτ_ref) < 1.0e-2
end

@testset "Few-mode bath, mixed time: all nine Green's functions" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = collect(0:δt:(Nt*δt))
	τs = collect(0:δτ:β)
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)

	# ED references (U = 0: the single-particle level is exact). The nine
	# branch pairs (b1, b2) ∈ {τ,+,−}² map onto the four elementary real-time
	# correlators plus the Matsubara one, with sign conventions calibrated
	# against GTEMPO's raw gf (annihilator-first conj on τ/+, creator-first on −):
	#   (τ,·)          →  ⟨d(τ_j)d†(0)⟩            = gτ_ed[j+1]
	#   (+,+)          →  ⟨d(t_j)d†(0)⟩            = R1[j]
	#   (−,+)          →  ⟨d†(0)d(t_j)⟩            = R2[j]
	#   (+,−)          →  −⟨d†(0)d(t_j)⟩           = −R2[j]
	#   (−,−)          →  −⟨d(0)d†(t_j)⟩           = −R3[j]
	#   (+,τ)          →  −⟨d†(0)d(t_j)⟩           = −R2[j]
	#   (−,τ)          →  −⟨d(0)d†(t_j)⟩           = −R3[j]
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	ρ = exp(-β * H)
	ρ /= tr(ρ)
	cache_ed = eigencache(H)
	R1 = correlation_2op_1t(H, a, adag, ρ, ts, cache_ed, reverse=false)
	R2 = correlation_2op_1t(H, adag, a, ρ, ts, cache_ed, reverse=true)
	R3 = correlation_2op_1t(H, a, adag, ρ, ts, cache_ed, reverse=true)
	gτ_ref = gτ_ed(H, a, adag, τs, β)

	spec = DiracDelta(ω=ω, α=α)
	bath = fermionicbath(spec, β=β)
	model = AndersonIM(U=0, μ=μ)
	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	cache = environments(lat, mpsK, mpsI)

	# raw GTEMPO correlator for a branch pair
	function rawgf(b1, i1, b2, i2)
		cached_gf(lat, (ContourIndex(i1, conj=(b1 == :-), branch=b1, band=1),
		                ContourIndex(i2, conj=!(b1 == :-), branch=b2, band=1)), mpsK, mpsI; cache=cache)
	end

	tol = 3.0e-2
	# τ as first branch: G(τ_j, 0) against the Matsubara reference (j = 1..Nτ,
	# i.e. τ = δτ..β−δτ; the KMS endpoints j = 0, Nτ+1 are covered by gτ below)
	for j in 1:Nτ
		@test _error(rawgf(:τ, j, :τ, 1), gτ_ref[j+1], tol) < 1
		@test _error(rawgf(:τ, j, :+, 1), gτ_ref[j+1], tol) < 1
		@test _error(rawgf(:τ, j, :-, 1), gτ_ref[j+1], tol) < 1
	end
	# real branches against the four elementary correlators
	for j in 1:Nt+1
		@test _error(rawgf(:+, j, :+, 1), R1[j], tol) < 1
		@test _error(rawgf(:-, j, :+, 1), R2[j], tol) < 1
		@test _error(rawgf(:+, j, :-, 1), -R2[j], tol) < 1
		@test _error(rawgf(:-, j, :-, 1), -R3[j], tol) < 1
		@test _error(rawgf(:+, j, :τ, 1), -R2[j], tol) < 1
		@test _error(rawgf(:-, j, :τ, 1), -R3[j], tol) < 1
	end
end
