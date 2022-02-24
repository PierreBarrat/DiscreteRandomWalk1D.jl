module DiscreteRandomWalk1D

using Distributions
using Parameters
using StatsBase

export RW
export step!

const RW = DiscreteRandomWalk1D

include("objects.jl")
include("walk.jl")

end # module
