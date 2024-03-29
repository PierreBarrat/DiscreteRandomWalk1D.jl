"""
	step!(rw::RandomWalker)

Advance `rw` by one step.
"""
function step!(rw::RandomWalker)
	@error "No `step!` method defined for $(typeof(rw))."
end

"""
	trajectory!(rw::RandomWalker, N::Int, ABC::Vararg{AbsorbingBC}; Δn=1)

Return trajectory of random walker for `N` steps, sampled every `Δn`.
Stop if the walk meets one of the absorbing boundary conditions `ABC`.
The return value is a tuple containing
1. the vector of positions of the walker
2. the corresponding time steps
3. the absorbing boundary that ended the walk if any, nothing otherwise
"""
function trajectory!(rw::RandomWalker, N::Int, ABCs::Vararg{AbsorbingBC}; Δn=1)
	X = Vector{Float64}(undef, div(N, Δn) + 1)
	local abc = nothing
	n = 0
	i = 1
	while n <= N
		X[i] = position(rw)
		i += 1
		# Check whether `rw` found an absorbing condition
		f, abc = isabsorbed(rw, ABCs...)
		f && break
		# Step
		for j in 1:Δn
			step!(rw)
		end
		n += Δn
	end

	return (x = X[1:(i-1)], t = 0:Δn:min(N,n), final = abc)
end
function trajectory!(rw::RandomWalker, N::Int; Δn=1)
	X = Vector{Float64}(undef, div(N, Δn) + 1)
	n = 0
	i = 1
	while n <= N
		X[i] = position(rw)
		for j in 1:Δn
			step!(rw)
		end
		i += 1
		n += Δn
	end

	return (x=X, t=0:Δn:N, final=nothing)
end
