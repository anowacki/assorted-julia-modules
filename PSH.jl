"""
# PSH

Module PSH converts seismograms in the ZNE frame into the PSH frame using the free-surface
transform.
"""
module PSH

using SAC

export psh, psh!

"""
    psh_matrix!(α, β, p, M) -> M

Fill in the PSH transfer matrix `M`, given the P-wave velocity `α`, S-wave velocity `β`
and the slowness `p`.  Velocities are in km/s, slowness in s/°.
"""
function psh_matrix!{T}(a, b, p, M::AbstractArray{T,2})
    size(M) == (3,3) || error("Transfer matrix must have size (3,3)")
    p /= 111.1 # Convert to s/km
    isevanescent(a, b, p) &&
        warn("psh_matrix!: Evanescent wave (vp=$a, vs=$b, p=$p)")
    qa = sqrt(1/a^2 - p^2)
    qb = sqrt(1/b^2 - p^2)
    # Kennett, GJI, 1991
    Vpz = -(1 - 2b^2*p^2)/(2a*qa)
    Vpr = p*b^2/a
    Vsz = p*b
    Vsr = (1 - 2b^2*p^2)/(2b*qb)
    Vht = 1/2
    M[:,:] .= [Vpz Vsz  0 
               Vpr Vsr  0 
                0   0  Vht]
end

"""
    psh_matrix(α, β, p) -> M

Return the PSH transfer matrix M, given the P-wave velocity `α`, S-wave velocity `β`
and the slowness `p`.  Velocities are in km/s, slowness in s/°
"""
psh_matrix(a, b, p) = psh_matrix!(a, b, p, Array{Complex{Float64}}(3,3))


"""
    isevanescent(α, β, p) -> ::Bool

Return true if the combination of Vp `α`, Vs `β` and slowness `p` lead to an evanescent
wave at the surface.
"""
function isevanescent(a, b, p)
    p = Float64(p/111.)
    a, b = Float64(a), Float64(b)
    1/a^2 - p^2 <= 0.0 || 1/b^2 - p^2 <= 0.0
end


"""
    psh!(Z, R, T, α, β, p) -> P, S, H

Converts three arrays, `SACtr`s or arrays of `SACtr`s in ZRT orientation into
PSH orientation, using the supplied Vp `α`, Vs `β` (km/s) and slowness `p` (s/°).

Traces are updated to be in the order P (like vertical), S (like radial), H (like transverse).
"""
function psh!{F}(Z::Array{F,1}, R::Array{F,1}, T::Array{F,1}, a, b, p)
	length(Z) == length(R) == length(T) || error("psh!: Traces must be same length")
	M = psh_matrix(a, b, p)
    # Kennett 1991
    Z, R, T = M[1,1]*Z + M[1,2]*R + M[1,3]*T,
              M[2,1]*Z + M[2,2]*R + M[2,3]*T,
              M[3,1]*Z + M[3,2]*R + M[3,3]*T
    Z, R, T
end
function psh!(Z::SACtr, R::SACtr, T::SACtr, a, b, p)
    psh!(Z.t, R.t, T.t, a, b, p)
    for t in (Z, R, T) SAC.update_headers!(t) end
    Z.kcmpnm, R.kcmpnm, T.kcmpnm = "P", "S", "H"
    Z, R, T
end
function psh!(Z::Array{SACtr}, R::Array{SACtr}, T::Array{SACtr}, a, b, p)
    for (z, r, t) in zip(Z, R, T) psh!(z, r, t, a, b, p) end
    Z, R, T
end

"""
    psh(Z, R, T, α, β, p) -> P, S, H

Convert three arrays, `SACtr`s or arrays of `SACtr`s in ZRT orientation into
PSH orientation, using the supplied Vp `α`, Vs `β` (km/s) and slowness `p` (s/°),
where `P` is parallel to P-wave motion positive away from the source, `S` is parallel
to SV and `H` is the transverse direction.
"""
psh(Z, R, T, a, b, p) = psh!(deepcopy(Z), deepcopy(R), deepcopy(T), a, b, p)

end # module
