### Fig.3 of "Efficient method for quantum impurity problems out of equilibrium"
const mpath = dirname(dirname(dirname(dirname(Base.@__DIR__))))

# println(mpath)

push!(LOAD_PATH, mpath * "/InfiniteDMRG/src")
push!(LOAD_PATH, mpath * "/GeneralHamiltonians/src")
push!(LOAD_PATH, mpath * "/TEBD/src")
push!(LOAD_PATH, mpath * "/FermionicTCMPS/src")

using FermionicTCMPS
using DelimitedFiles, JSON

function J(D::Real, ω::Real)
	t′ = 0.3162 
	# return (sqrt(4*t^2-ω^2) / (2*π*t^2)) * t1^2 
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * t′^2 / 2
	# return (2/(D*pi)) * sqrt(1 - (ω/D)^2 ) * t′^2 
end


spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

function main_ed(V_over_Gamma, tall=6.3, dw=0.001)
	println("run for V = ", V_over_Gamma, " dw=", dw)
	β = Inf
	D = 2

	Γ = 0.1
	V = V_over_Gamma * Γ

	leftmu = V / 2
	rightmu = -V / 2
	ϵ_d = 0

	δt = tall / Γ
	t = tall / Γ
	# t = 0.7 / Γ
	N = round(Int, t / δt)
	
	

	leftbath = fermionicbath(spectrum_func(D), β=β, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(D), β=β, μ=rightmu)
	leftconfig = star(leftbath, dw=dw)
	rightconfig = star(rightbath, dw=dw)

	model = FreeSIDBD(leftconfig, rightconfig, μ = ϵ_d)
	ham = FermionicCommutator(coefficient_matrix(FreeFermionicHamiltonian(model)))
	observer1 = coefficient_matrix(FreeFermionicHamiltonian(left_sysbath_tunneling(model)), length(model))
	observer2 = coefficient_matrix(FreeFermionicHamiltonian(right_sysbath_tunneling(model)), length(model))

	ρ₀ = separable_state(model, sys_states=[0])
	sys_site = only(default_sys_sites(model))
	ns = [real(ρ₀[sys_site, sys_site])]
	currents1 = [sum(observer1 .* ρ₀)]
	currents2 = [sum(observer2 .* ρ₀)]
	for i in 1:N
		stepper = ExactStepper(tspan=(-im*(i-1)*δt, -im*i*δt), ishermitian=true)
		timeevo!(ρ₀, ham, stepper)
		push!(currents1, sum(observer1 .* ρ₀))
		push!(currents2, sum(observer2 .* ρ₀))
		push!(ns, real(ρ₀[sys_site, sys_site]))
	end	
	currents1 = -2*im .* currents1
	currents2 = -2*im .* currents2

	ts = [i*δt for i in 0:N]

	# return  ns, real(currents1), real(currents2)

	data_path = "result/anderson_ed_V$(V_over_Gamma)_t$(tall)_dw$(dw).json"

	results = Dict("ts"=>ts, "ns"=>ns, "Ileft" => real(currents1), "Iright"=>real(currents2))

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end
end

function main_ed_all(tall=6.3, dw=0.001)
	Vs = [0.17857143, 0.35714286, 0.53571429, 0.71428571, 0.89285715]
	for V in Vs
		main_ed(V, tall, dw)
	end
end