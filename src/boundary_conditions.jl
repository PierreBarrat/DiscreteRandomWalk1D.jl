abstract type BoundaryCondition end

######################################################################
######################### Boundary conditions ########################
######################################################################
const ABOVE = (:above, :right, :high, :up)
const BELOW = (:below, :left, :low, :down)

"""
	struct AbsorbingBC{T}

Absorbing boundary condition for a `RandomWalker`. The field `sense` can be \
`$(ABOVE)` or `$(BELOW)`.

A boundary condition `bc` works like this:
- if `sense` is `:above` or similar, then `position(walker) >= bc.A` will trigger absorption
- if `sense` is `:below` or similar, then `position(walker) <= bc.A` will trigger absorption

## Fields

```
A::T
sense::Symbol
```
"""
@Base.kwdef struct AbsorbingBC{T} <: BoundaryCondition
	A::T
	sense::Symbol
end


function isabsorbed(rw::RandomWalker, ABCs::Vararg{AbsorbingBC})
	for abc in ABCs
		if isabsorbed(rw, abc)[1]
			return true, abc
		end
	end
	return false, nothing
end
function isabsorbed(rw::RandomWalker, ABC::AbsorbingBC)
	if in(ABC.sense, ABOVE) && position(rw) >= ABC.A
		return true, ABC
	elseif in(ABC.sense, BELOW) && position(rw) <= ABC.A
		return true, ABC
	else
		return false, nothing
	end
end

natural_boundaries(::NeutralCoalescent) = AbsorbingBC(1, :down)
natural_boundaries(::EFCoalescent) = AbsorbingBC(1, :down)
