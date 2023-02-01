### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 88e66206-a21e-11ed-119f-531ed9b4a18f
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(; temp=true)
	Pkg.add(url="https://github.com/PierreBarrat/DiscreteRandomWalk1D.jl")
	Pkg.add("Plots")
	using DiscreteRandomWalk1D
	using Plots
end

# ╔═╡ 720232f0-39d9-498f-ad20-1e1148fc6c93
md"""
# Defining a random walker: `DrunkMan`

The drunk man walks in a street, starting at some initial position `x0`. He takes a step right with probability `p_right`, and a step left with probability `p_left`. The size of the steps is `dx = 0.5`. meters. 
"""

# ╔═╡ fd17e08f-27f2-43ae-97bd-b7c6c4baa902
"""
	mutable struct DrunkMan

Type representing the drunk man
"""
Base.@kwdef mutable struct DrunkMan <: RandomWalker
	p_right::Float64 = 0.5
	dx::Float64 = 0.5
	t::Int = 0
	x::Float64 = 0.
end

# ╔═╡ aa314e5f-757a-4fec-b753-9adcb0bf87f7
md"""
We now need to extend  the methods `step!` and `position` from `DiscreteRandomWalk1D`:
- `step!(rw::RandomWalker)` makes `rw` take one step. The output value is ignored.
- `position(rw::RandomWalker` returns the scalar position. 
"""

# ╔═╡ bc8caf69-ac08-4923-9711-b8c481944e14
begin
	import DiscreteRandomWalk1D: step!, position
	position(rw::DrunkMan) = rw.x
	
	function step!(rw::DrunkMan)
		if rand() < rw.p_right
			rw.x += rw.dx
		else
			rw.x -= rw.dx
		end
		rw.t += 1
	
		return nothing
	end
end

# ╔═╡ ac1c98bd-cddc-4683-ba9b-f9e588e15095
md"""
Now, let's have the drunk man take a walk. For this we use the `trajectory!` function. We then plot the trajectory as a function of time. 
The following parameters are used: 
- initial position `x0 = $(x0)`
- probability of a step to the right: `p_right = $(p_right)`
"""

# ╔═╡ a1319200-aa76-40a5-8a9f-754582c532d9
begin
	x0 = 10
	p_right = 0.52
end;

# ╔═╡ a77d052c-04de-47bd-bcb7-8dc16e68c493
X = let
	drunkman = DrunkMan(; x = x0, p_right)
	RW.trajectory!(drunkman, 1000)
end

# ╔═╡ 213c7d2f-b0a6-440d-b6e4-aa8effa96e20
md"""
`X` has three fields: 
- `X.x` is a vector of positions
- `X.t` is a vector of time points
- `X.final` handles the final state of the random walker. It is useful when using boundary conditions. 
"""

# ╔═╡ 85c385c3-c98a-4098-acc0-c2f5e111a8bd
let
	p = plot(
		X.t, X.x;
		label="", title = "Position of the drunkman", xlabel="time", ylabel="x"
	)
end

# ╔═╡ d78a6261-63b8-4f72-8512-6ac987235041
md"""
The drunkman is trying to reach the bar, at position `x=20`: this is why he is trying to step right (`p_right > 0.5`). However, what he really needs is to get home, at position `x=0`. 

To represent this, we use absorbing boundary conditions.
"""

# ╔═╡ 9a6c4b55-318c-496d-ad4c-1ff8ce454730
begin
	home = RW.AbsorbingBC(A = 0, sense = :left)
	bar = RW.AbsorbingBC(A = 20, sense = :right)
end;

# ╔═╡ 510454c1-12cc-4b61-b673-5265246448d8
md"""
If the drunk man ever goes left (or below, depending on how you see it) `home.A`, he will be "absorbed" by  the boundary condition: he made it home!
If he goes right of `bar.A`, he reached the bar. 

Boundary conditions are simply added as extra arguments to `trajectory!`. 
The output `X` of trajectory stores the final state in the field `X.final`.
"""

# ╔═╡ ee379fbe-e50a-44c9-b1b9-dd20d1388b19
X_with_boundaries = let
	drunkman = DrunkMan(; x = x0, p_right)
	RW.trajectory!(drunkman, 1_000, home, bar)
end

# ╔═╡ a2270eaf-d580-4d4d-8c7f-79ba9e261113
if X_with_boundaries.final == home
	println("Made it home! Hourra!")	
elseif X_with_boundaries.final == bar
	println("Reached the bar! Time for another drink!")	
elseif isnothing(X_with_boundaries.final)
	println("Still wandering in the streets...")
end

# ╔═╡ e785f3bc-b2c2-44da-a08c-8ea24decdd4d
let
	p = plot(
		X_with_boundaries.t, X_with_boundaries.x;
		label="", title = "Position of the drunkman", xlabel="time", ylabel="x"
	)
	hline!([home.A], line=(3, :green, :dash), label="Home")
	hline!([bar.A], line=(3, :red, :dash), label="Bar")
	plot!(legend=:outerright)
end

# ╔═╡ Cell order:
# ╠═88e66206-a21e-11ed-119f-531ed9b4a18f
# ╟─720232f0-39d9-498f-ad20-1e1148fc6c93
# ╠═fd17e08f-27f2-43ae-97bd-b7c6c4baa902
# ╟─aa314e5f-757a-4fec-b753-9adcb0bf87f7
# ╠═bc8caf69-ac08-4923-9711-b8c481944e14
# ╟─ac1c98bd-cddc-4683-ba9b-f9e588e15095
# ╠═a1319200-aa76-40a5-8a9f-754582c532d9
# ╠═a77d052c-04de-47bd-bcb7-8dc16e68c493
# ╟─213c7d2f-b0a6-440d-b6e4-aa8effa96e20
# ╟─85c385c3-c98a-4098-acc0-c2f5e111a8bd
# ╟─d78a6261-63b8-4f72-8512-6ac987235041
# ╠═9a6c4b55-318c-496d-ad4c-1ff8ce454730
# ╟─510454c1-12cc-4b61-b673-5265246448d8
# ╠═ee379fbe-e50a-44c9-b1b9-dd20d1388b19
# ╠═a2270eaf-d580-4d4d-8c7f-79ba9e261113
# ╟─e785f3bc-b2c2-44da-a08c-8ea24decdd4d
