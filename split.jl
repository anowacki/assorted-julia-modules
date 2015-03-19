module split
# Perform shear wave splitting analysis on two orthogonal SAC traces

using SAC
SACtr = SAC.SACtr

# Module defaults which can be changed
split_ndt = 51
split_nphi = 181
split_dt_max = 4.

function search_phi(t1::SACtr, t2::SACtr;
					nphi=split_nphi, ndt=split_ndt, dt_max=split_dt_max)
	# Return the array of minimum eigenvalues of the covariance matrices
	# when examining trial operators with nphi fast orientations
	phi = linspace(-90., 90., nphi)
	dt = zeros(ndt)
	lam1, lam2 = zeros(nphi,ndt), zeros(nphi,ndt)
	for i = 1:nphi
		dt, lam1[i,:], lam2[i,:] = search_dt(t1, t2, phi[i], dt_max, ndt)
	end
	return phi, dt, lam1, lam2
end

function search_dt(t1::SACtr, t2::SACtr, phi, dt_max, ndt)
	# For a given fast orientation, return the smaller eigenvalue of the
	# covariance matrix in an array accompanied by a vector of dts.
	dt = linspace(dt_max/ndt, dt_max, ndt)
	lam1, lam2 = zeros(ndt), zeros(ndt)
	for i = 1:ndt
		lam1[i], lam2[i] = compute_eigvals(t1, t2, phi, dt[i])
	end
	return dt, lam1, lam2
end

function apply_split!(t1::SACtr, t2::SACtr, phi::Number, dt::Number)
	# Apply a splitting operator to two SAC traces.  They are assumed to be
	# orthogonal, and the first trace should be anticlockwise of the second.
	theta = phi - t1.cmpaz
	SAC.rotate_through!(t1, t2, theta)
	# Move the fast trace forward
	SAC.tshift!(t1, -dt)
	SAC.rotate_through!(t1, t2, -theta)
	return
end

function compute_eigvals(t1::SACtr, t2::SACtr, phi::Number, dt::Number)
	# Compute the two eigenvalues for traces for a specific fast orientation, phi,
	# and delay time, dt.
	T1, T2 = SAC.copy(t1), SAC.copy(t2)
	apply_split!(T1, T2, phi, -dt)
	c = covar(T1, T2)
	eval = Base.eigvals(c)
	return maximum(eval), minimum(eval)
end

function covar(t1::SACtr, t2::SACtr)
	# Return the covariance matrix for two traces
	c = Array{SAC.SACFloat}(2,2)
	c[1,1] = sum(t1.t.^2)
	c[2,2] = sum(t2.t.^2)
	c[1,2] = c[2,1] = sum(t1.t.*t2.t)
	return c
end

end # module
