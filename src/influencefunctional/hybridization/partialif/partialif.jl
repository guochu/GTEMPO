include("util.jl")

# static mode, compute the whole IF and then measure
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")
include("partialif_static.jl")

# dynamic mode, evolve the IF and measure during evolve, similar to TEMPO
include("realtime_stepper.jl")
include("partialif_stepper.jl")