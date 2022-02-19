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
γ::Float64 = exp(-s/α)
t::Float64 = 0.
x::Float64 = 0.
```
"""
@with_kw mutable struct EFRW <: RandomWalker
	s::Float64 = 0.01
	α::Float64 = 0.01
	γ::Float64 = exp(-s/α)
	t::Float64 = 0.
	x::Float64 = 0.
end

p_right(rw::EFRW) = 1/2
p_left(rw::EFRW) = 1/2
step_size_right(rw::EFRW) = (1-rw.γ) * (1-rw.x)
step_size_left(rw::EFRW) = (1-rw.γ) * rw.x



"""
	mutable struct EFRW

Approximation to `EFRW` with a symetric step size.

## Fields

```
s::Float64 = 0.01
α::Float64 = 0.01
γ::Float64 = exp(-s/α)
t::Float64 = 0.
x::Float64 = 0.
```
"""
@with_kw mutable struct SymEFRW <: RandomWalker
	s::Float64 = 0.01
	α::Float64 = 0.01
	γ::Float64 = exp(-s/α)
	#
	t::Float64 = 0.
	x::Float64 = 0.
end

p_right(rw::SymEFRW) = (1 + (1-2*rw.x)/sqrt(2*((1-rw.x)^2 + rw.x^2)))/2
p_left(rw::SymEFRW) = 1. - p_right(rw)

step_size_right(rw::SymEFRW) = (1-rw.γ) * sqrt(((1-rw.x)^2 + rw.x^2)/2)
step_size_left(rw::SymEFRW) = step_size_right(rw)
