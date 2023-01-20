"""
	step!(rw::RandomWalker)

Advance `rw` by one step.
Call `p_right` (resp. left) and `step_size_right` (resp. left).
"""
function step!(rw::RandomWalker)
	@error "No `step!` method defined for $(typeof(rw))."
end
# function step!(rw::RandomWalker)
# 	p = rand()
# 	if p < p_right(rw)
# 		rw.x += step_size_right(rw)
# 	else
# 		rw.x -= step_size_left(rw)
# 	end
# 	update!(rw)
# 	rw.t += 1
# 	return rw
# end

"""
	trajectory!(rw, N, ABC::Vararg{AbsorbingBC}; Δn=1)

Return trajectory of random walker for `N` steps, sampled every `Δn`.
Stop if the walk meets one of the absorbing boundary conditions `ABC`.
The return value is a tuple containing
1. the vector of positions of the walker
2. the corresponding time steps
3. the absorbing boundary that ended the walk if any, nothing otherwise
"""
function trajectory!(rw, N, ABC, ABCs::Vararg{AbsorbingBC}; Δn=1)
	X = Vector{Float64}(undef, Int(floor(N/Δn)) + 1)
	local abc = nothing
	n = 0
	i = 1
	while n <= N
		X[i] = position(rw)
		i += 1
		# Check whether `rw` found an absorbing condition
		f, abc = isabsorbed(rw, ABC, ABCs...)
		f && break
		# Step
		for j in 1:Δn
			step!(rw)
		end
		n += Δn
	end

	return (x = X[1:(i-1)], t = 0:Δn:min(N,n), final = abc)
end
function trajectory!(rw, N; Δn=1)
	X = Vector{Float64}(undef, Int(floor(N/Δn)) + 1)
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

	return X, 0:Δn:N
end
