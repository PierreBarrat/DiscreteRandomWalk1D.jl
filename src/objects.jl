abstract type RandomWalker end

######################################################################
########################### EF random walks ##########################
######################################################################

"""
	mutable struct EFRW

Expiring Fitness Random Walker.

## Fields

```
β::Distribution = Dirac(0.5)
ρ::Float64 = 1.
μ::Float64 = 0.
N::Float64 = 0
t::Int = 0
x::Float64 = 0.
```
"""
Base.@kwdef mutable struct EFRW <: RandomWalker
	β::Distribution = Dirac(0.5)
	ρ::Float64 = 1.
	μ::Float64 = 0.
	N::Float64 = 0
	t::Int = 0
	x::Float64 = 0.
end
"""
	EFRW(β; kwargs...)

Create an `EFRW` object with the amplitude of partial sweeps distributed according to `β`.
If `β::Number`, the distribution will be `Dirac(β)`.

For a description of keyword arguments, see the help for the struct `EFRW`.
"""
EFRW(β::Distribution; kwargs...) = EFRW(; β, kwargs...)
EFRW(β::Number; kwargs...) = EFRW(; β = Dirac(β), kwargs...)

function step!(rw::EFRW)
	if rw.μ > 0
		if rw.N > 0
			nx = round(rw.N * rw.x) # number of individuals carrying x
			tx = rand(Binomial(rw.N-nx, rw.μ))/rw.N # mutations to x
			fx = rand(Binomial(nx, rw.μ))/rw.N # mutations from x
			rw.x += tx - fx
		else
			# N infinite: mutations are deterministic
			rw.x += rw.μ * (1 - 2*rw.x)
		end
	end

	# Partial sweep
	if rw.ρ >= 1 || rand() < rw.ρ # Is there a partial sweep?
		if rand() <= rw.x
			rw.x += rand(rw.β) * (1. - rw.x)
		else
			rw.x -= rand(rw.β) * rw.x
		end
	end

	# Binomial sampling if finite population size
	if rw.N > 0
		rw.x = rand(Binomial(rw.N, rw.x)) / rw.N
	end

	rw.t += 1

	return nothing
end

position(rw::EFRW) = rw.x
pos(rw::EFRW) = rw.x


######################################################################
######################### Neutral random walk ########################
######################################################################

Base.@kwdef mutable struct NeutralRW <: RandomWalker
	N::Int = 100
	t::Int = 0
	x::Float64 = 0.5
end

function step!(rw::NeutralRW)
	rw.x = rand(Binomial(rw.N, rw.x)) / rw.N
	rw.t += 1
	return nothing
end

position(rw::NeutralRW) = rw.x
pos(rw::NeutralRW) = rw.x

######################################################################
############################# Coalescent #############################
######################################################################

"""
	NeutralCoalescent(;
		N::Int = 1000
		t::Int = 0
		n::Int = 10
	)
"""
Base.@kwdef mutable struct NeutralCoalescent <: RandomWalker
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

	return nothing
end

position(rw::NeutralCoalescent) = rw.n
pos(rw::NeutralCoalescent) = position(rw)


"""
	EFCoalescent(;
		N::Float64 = Inf
		β::Float64 = 0.5
		ρ::Float64 = 0.05
		t::Int = 0
		n::Int = 10
	)
"""
Base.@kwdef mutable struct EFCoalescent <: RandomWalker
	N::Float64 = Inf
	β::Float64 = 0.5
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

	return nothing
end

position(rw::EFCoalescent) = rw.n
pos(rw::EFCoalescent) = position(rw)


