######################################################################
######################### Boundary conditions ########################
######################################################################
const ABOVE = (:above, :right, :high, :up)
const BELOW = (:below, :left, :low, :down)

"""
	struct AbsorbingBC{T}

Absorbing boundary condition for a `RandomWalker`. The field `sense` can be \
`$(ABOVE)` or `$(BELOW)`.

## Fields

```
A::T
sense::Symbol
```
"""
struct AbsorbingBC{T}
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
