const DefaultMPOTruncation = truncdimcutoff(D=10000, ϵ=1.0e-10, add_back=0)

# function grassmanncreation()
# 	t = TensorMap(ds->zeros(Float64, ds), _ph ⊗ _ph, _ph ⊗ _ph)
# 	f1 = FusionTree((Z2Irrep(1), Z2Irrep(1)), Z2Irrep(0), (false, false))
# 	f2 = FusionTree((Z2Irrep(0), Z2Irrep(0)), Z2Irrep(0), (false, false))
# 	t[f1, f2] .= 1	
# 	return t
# end

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
