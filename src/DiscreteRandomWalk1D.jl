"""
	module DiscreteRandomWalk1D

## Example

```jldoctest
julia> walker = RW.NeutralRW(; N = 100, x = 0.5) # neutral evolution, starting at x=0.5
DiscreteRandomWalk1D.NeutralRW
  N: Int64 100
  t: Int64 0
  x: Float64 0.5


julia> boundaries = [RW.AbsorbingBC(0.05, :left), RW.AbsorbingBC(0.95, :right)]
2-element Vector{DiscreteRandomWalk1D.AbsorbingBC{Float64}}:

 DiscreteRandomWalk1D.AbsorbingBC{Float64}(0.05, :left)
 DiscreteRandomWalk1D.AbsorbingBC{Float64}(0.95, :right)

julia> X = RW.trajectory!(
               walker,
               1000, # 1000 time steps
               boundaries...;
               Δn = 10 # record position every 10th step only
       )
([0.5, 0.4, 0.19, 0.5, 0.7, 0.88, 0.91, 1.0], 0:10:70, DiscreteRandomWalk1D.AbsorbingBC{Float64}(0.95, :right))

```

## Walkers

Subtypes of `RandomWalker`:

- `NeutralRW`: Neutral evolution with binomial sampling
- `EFRW`: expiring fitness
- `NeutralCoalescent`: number of lineages in a Kingman coalescent
- `EFCoalescent`: number of lineages in an expiring fitness coalescent

## Exported names

RW, step!
"""
module DiscreteRandomWalk1D

using Distributions

export RW

export RandomWalker
export AbsorbingBC

export step!, trajectory

const RW = DiscreteRandomWalk1D

include("objects.jl")
include("boundary_conditions.jl")
include("walk.jl")

end # module
