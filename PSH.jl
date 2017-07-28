"""
# PSH

Module PSH converts seismograms in the ZNE frame into the PSH frame using the free-surface
transform.
"""
module PSH

using SAC

export psh!

"""
    psh_matrix!(α, β, p, M)

Fill in the PSH transfer matrix `M`, given the P-wave velocity `α`, S-wave velocity `β`
and the slowness `p`.  Velocities are in km/s, slowness in s/°
"""
function psh_matrix!{T}(a, b, p, M::AbstractArray{T,2})
    size(M) == (3,3) || error("Transfer matrix must have size (3,3)")
    pr = p/111.
    isevanescent(a, b, p) &&
		error("psh_matrix!: Evanescent wave not permitted (vp=$a, vs=$b, p=$p)")
    qa = sqrt(1/a^2 - pr^2)
    qb = sqrt(1/b^2 - pr^2)
    M[:,:] .= [pr*b^2/a              0   (b^2*pr^2-0.5)/(a*qa);
               (0.5-b^2*pr^2)/(b*qb) 0   pr*b;
               0                    0.5  0                    ]
end

"""
    psh_matrix(α, β, p) -> M

Return the PSH transfer matrix M, given the P-wave velocity `α`, S-wave velocity `β`
and the slowness `p`.  Velocities are in km/s, slowness in s/°
"""
function psh_matrix(a, b, p)
    M = Array{Float64}(3, 3)
    psh_matrix!(a, b, p, M)
    M
end

"""
    isevanescent(α, β, p) -> ::Bool

Return true if the combination of Vp `α`, Vs `β` and slowness `p` lead to an evanescent
wave at the surface.
"""
isevanescent(a, b, p) = 1/a^2 - (p/111.)^2 <= 0 || 1/b^2 - (p/111.)^2 <= 0


"""
    psh!(Z, R, T, α, β, p) -> P, S, H

Converts three arrays, `SACtr`s or arrays of `SACtr`s in ZRT orientation into
PSH orientation, using the supplied Vp `α`, Vs `β` (km/s) and slowness `p` (s/km).

Traces are updated to be in the order P (like vertical), S (like radial), H (like transverse).
"""
function psh!{F}(Z::Array{F,1}, R::Array{F,1}, T::Array{F,1}, a, b, p)
	length(Z) == length(R) == length(T) || error("psh!: Traces must be same length")
	M = psh_matrix(a, b, p)
    A = M*[R'; T'; Z']
    Z[:] = vec(A[1,:])
    R[:] = vec(A[2,:])
    T[:] = vec(A[3,:])
    Z, R, T
end
function psh!(Z::SACtr, R::SACtr, T::SACtr, a, b, p)
    psh!(Z.t, R.t, T.t, a, b, p)
    for t in (Z, R, T) SAC.update_headers!(t) end
    Z, R, T
end
function psh!(Z::Array{SACtr}, R::Array{SACtr}, T::Array{SACtr}, a, b, p)
    for (z, r, t) in zip(Z, R, T) psh!(z, r, t) end
    Z, R, T
end

"""
    psh(Z, R, T, α, β, p) -> P, S, H

Convert three arrays, `SACtr`s or arrays of `SACtr`s in ZRT orientation into
PSH orientation, using the supplied Vp `α`, Vs `β` (km/s) and slowness `p` (s/km),
where `P` is parallel to P-wave motion positive away from the source, `S` is parallel
to SV and `H` is the transverse direction.
"""
psh(Z, R, T, a, b, p) = psh!(deepcopy(Z), deepcopy(R), deepcopy(T), a, b, p)

end # module
