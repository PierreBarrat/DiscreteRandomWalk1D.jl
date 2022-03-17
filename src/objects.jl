abstract type RandomWalker end

######################################################################
########################### EF random walks ##########################
######################################################################

"""
	mutable struct EFRW

Expiring Fitness Random Walker.

## Fields

```
s::Float64 = 0.01
α::Float64 = 0.01
β::Float64 = 1 - exp(-s/α)
t::Float64 = 0.
x::Float64 = 0.
```
"""
@with_kw mutable struct EFRW <: RandomWalker
	s::Float64 = 0.01
	α::Float64 = 0.01
	β::Float64 = 1. - exp(-s/α)
	t::Int = 0
	x::Float64 = 0.
end

p_right(rw::EFRW) = rw.x
p_left(rw::EFRW) = 1. - rw.x
step_size_right(rw::EFRW) = rw.β * (1. - rw.x)
step_size_left(rw::EFRW) = rw.β * rw.x

function step!(rw::EFRW)
	p = rand()
	if p <= rw.x
		rw.x += rw.β * (1. - rw.x)
	else
		rw.x -= rw.β * rw.x
	end
	rw.t += 1
	return rw
end

position(rw::EFRW) = rw.x
pos(rw::EFRW) = rw.x

"""
	mutable struct StochasticEFRW

Expiring Fitness Random Walker with a randomly chosen fitness effect at each step.

## Fields

```
s::Distribution
α::Float64 = 0.01
β::Float64 = 1 - exp(-s/α)
t::Int = 0
x::Float64 = 0.
```
"""
@with_kw mutable struct StochasticEFRW <: RandomWalker
	s::Distribution
	sval::Float64 = rand(s)
	α::Float64 = 0.01
	β::Float64 = 1. - exp(-sval/α)
	t::Int = 0
	x::Float64 = 0.
end

p_right(rw::StochasticEFRW) = rw.x
p_left(rw::StochasticEFRW) = 1. - rw.x
step_size_right(rw::StochasticEFRW) = rw.β * (1. - rw.x)
step_size_left(rw::StochasticEFRW) = rw.β * rw.x

function step!(rw::StochasticEFRW)
	p = rand()
	if p <= rw.x
		rw.x += rw.β * (1. - rw.x)
	else
		rw.x -= rw.β * rw.x
	end
	rw.t += 1
	# Update β
	rw.sval = rand(rw.s)
	rw.β = 1. - exp(-rw.sval/rw.α)

	return rw
end

position(rw::StochasticEFRW) = rw.x
pos(rw::StochasticEFRW) = rw.x


# """
# 	mutable struct SymEFRW

# Approximation to `EFRW` with a symetric step size.

# ## Fields

# ```
# s::Float64 = 0.01
# α::Float64 = 0.01
# γ::Float64 = exp(-s/α)
# t::Int = 0
# x::Float64 = 0.
# ```
# """
# @with_kw mutable struct SymEFRW <: RandomWalker
# 	s::Float64 = 0.01
# 	α::Float64 = 0.01
# 	γ::Float64 = exp(-s/α)
# 	#
# 	t::Int = 0
# 	x::Float64 = 0.5
# end

# p_right(rw::SymEFRW) = (1 + (1-2*rw.x)/sqrt(2*((1-rw.x)^2 + rw.x^2)))/2
# p_left(rw::SymEFRW) = 1. - p_right(rw)

# step_size_right(rw::SymEFRW) = (1-rw.γ) * sqrt(((1-rw.x)^2 + rw.x^2)/2)
# step_size_left(rw::SymEFRW) = step_size_right(rw)

# update!(::SymEFRW) = nothing

######################################################################
######################### Neutral random walk ########################
######################################################################

@with_kw mutable struct NeutralRW <: RandomWalker
	N::Int = 100
	t::Int = 0
	x::Float64 = 0.5
end

p_right(::NeutralRW) = 1.
p_left(::NeutralRW) = 0.
function step_size_right(rw::NeutralRW)
	Y = Binomial(rw.N, rw.x)
	return rand(Y)/rw.N - rw.x
end
step_size_left(::NeutralRW) = 0.
update!(::NeutralRW) = nothing

function step!(rw::NeutralRW)
	rw.x = rand(Binomial(rw.N, rw.x)) / rw.N
	rw.t = 1
	return rw
end

position(rw::NeutralRW) = rw.x
pos(rw::NeutralRW) = rw.x

# ######################################################################
# ########################## Neutral genealogy #########################
# ######################################################################

# @with_kw mutable struct NeutralCoalescent <: RandomWalker
# 	N::Int = 1000
# 	t::Int = 0
# 	n::Int = 10
# end


######################################################################
######################### Boundary conditions ########################
######################################################################
const ABOVE = (:above, :right, :high, :up)
const BELOW = (:below, :left, :low, :down)

"""
	struct AbsorbingBC

Absorbing boundary condition for a `RandomWalker`. The field `sense` can be \
`$(ABOVE)` or `$(BELOW)`.

## Fields

```
A::Float64
sense::Symbol
```
"""
struct AbsorbingBC
	A::Float64
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
	if in(ABC.sense, ABOVE) && rw.x >= ABC.A
		return true, ABC
	elseif in(ABC.sense, BELOW) && rw.x <= ABC.A
		return true, ABC
	else
		return false, nothing
	end
end
