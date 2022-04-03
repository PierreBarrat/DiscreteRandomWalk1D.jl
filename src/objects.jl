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

Expiring Fitness Random Walker with a randomly chosen β.

## Fields

```
β::Distribution = (1. - exp(-sval/α)) * DiscreteUniform(1,1)
βval::Float64 = rand(β)
t::Int = 0
x::Float64 = 0.
```
"""
@with_kw mutable struct StochasticEFRW <: RandomWalker
	β::Distribution = (1. - exp(-sval/α)) * DiscreteUniform(1,1)
	βval::Float64 = rand(β)
	t::Int = 0
	x::Float64 = 0.
end
function StochasticEFRW(β::Float64; kwargs...)
	return StochasticEFRW(β = β*DiscreteUniform(1,1); kwargs...)
end

function step!(rw::StochasticEFRW)
	p = rand()
	if p <= rw.x
		rw.x += rw.βval * (1. - rw.x)
	else
		rw.x -= rw.βval * rw.x
	end
	rw.t += 1
	# Update β
	rw.βval = rand(rw.β)

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

function step!(rw::NeutralRW)
	rw.x = rand(Binomial(rw.N, rw.x)) / rw.N
	rw.t = 1
	return rw
end

position(rw::NeutralRW) = rw.x
pos(rw::NeutralRW) = rw.x

# ######################################################################
# ############################# Coalescent #############################
# ######################################################################

"""
	NeutralCoalescent(;
		N::Int = 1000
		t::Int = 0
		n::Int = 10
	)
"""
@with_kw mutable struct NeutralCoalescent <: RandomWalker
	N::Int = 1000
	t::Int = 0
	n::Int = 10
end

#=
The coalescence *rate* is n(n-1)/2N
=#
function step!(rw::NeutralCoalescent)
	R = rw.n * (rw.n - 1) / 2 / rw.N
	# For low rate
	if R < 0.05
		p = rand()
		if p < rw.n * (rw.n - 1) / 2 / rw.N
			rw.n -= 1
		end
	else
		δn = rand(Poisson(R))
		rw.n -= δn
	end
	rw.t += 1
	# Minimum of one lineage
	if rw.n < 1
		rw.n = 1
	end

	return rw
end

position(rw::NeutralCoalescent) = rw.n
pos(rw::NeutralCoalescent) = position(rw)


"""
	EFCoalescent(;
		N::Float64 = Inf
		s::Float64 = 0.01
		α::Float64 = 0.01
		β::Float64 = 1. - exp(-s/α)
		ρ::Float64 = 0.05
		t::Int = 0
		n::Int = 10
	)
"""
@with_kw mutable struct EFCoalescent <: RandomWalker
	N::Float64 = Inf
	s::Float64 = 0.01
	α::Float64 = 0.01
	β::Float64 = 1. - exp(-s/α)
	ρ::Float64 = 0.05
	t::Int = 0
	n::Int = 10
end

function step!(rw::EFCoalescent)
	p = rand()
	if p < rw.n * (rw.n - 1) / 2 / rw.N
		rw.n -= 1
	elseif p < rw.n * (rw.n - 1) / 2 / rw.N + rw.ρ
		# On average βn lineages coalesce into 1
		R = rand(Binomial(rw.n, rw.β))
		if R > 0
			rw.n += 1-R
		end
	end
	rw.t += 1
	# Minimum of one lineage
	if rw.n < 1
		rw.n = 1
	end

	return rw
end

position(rw::EFCoalescent) = rw.n
pos(rw::EFCoalescent) = position(rw)

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
