"""
	step!(rw::RandomWalker)

Advance `rw` by one step.
Call `p_right` (resp. left) and `step_size_right` (resp. left).
"""
function step!(rw::RandomWalker)
	p = rand()
	if p < p_right(rw)
		rw.x += step_size_right(rw)
	else
		rw.x -= step_size_left(rw)
	end
	rw.t += 1.
	return rw
end

"""
	trajectory!(rw, N, n=1)

Return trajectory of random walker for `N` steps, sampled every `n`.
"""
function trajectory!(rw, N, Δn=1)
	X = Vector{Float64}(undef, Int(floor(N/Δn)) + 1)
	X[1] = rw.x

	n = Δn
	i = 2
	while n <= N
		for j in 1:Δn
			step!(rw)
		end
		X[i] = rw.x
		i += 1
		n += Δn
	end

	return X, 0:Δn:N
end
