# the functions here are only applicable for specific orderings

# There is a "root" ordering which can be used to build K accurately and efficiently
# TODO: generalize accsysdynamics to BandLocalLayout


include("util.jl")
include("zooming.jl")

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")