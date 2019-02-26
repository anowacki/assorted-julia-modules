module Split
# Perform shear wave splitting analysis on two orthogonal SAC traces

using LinearAlgebra
using StaticArrays

using SAC

# Module defaults
const split_ndt = 51
const split_nphi = 181
const split_dt_max = 4.

"""
    search_phi(t1, t2; nphi, ndt, dt_max) -> results

Perform a search over a pair of SAC traces, `t1` and `t2`, for the smallest value of the
minimum eigenvalue of the covariance matrix between the traces, for a set of `nphi`×`ndt`
shear wave splitting operators.

`results` is a named tuple containing:

- `phi`: The set of fast shear wave orientations in °
- `dt`: The set of delays times in s
- `lam1`: The larger eigenvalues at each [phi,dt] point
- `lam2`: The smaller eigenvalues at each point
- `phi_best`: The best ϕ
- `dt_best`: The best δt
- `spol` and `dspol`: The source polarisation and an estimate of its uncertainty for the
  best-fitting ϕ and δt
"""
function search_phi(t1::SACtr, t2::SACtr;
                    nphi=split_nphi, ndt=split_ndt, dt_max=split_dt_max)
    # Return the array of minimum eigenvalues of the covariance matrices
    # when examining trial operators with nphi fast orientations
    phi = range(-90., stop=90., length=nphi)
    dt = range(dt_max/ndt, stop=dt_max, length=ndt)
    lam1, lam2 = zeros(SAC.SACFloat, nphi,ndt), zeros(SAC.SACFloat, nphi,ndt)
    T1, T2 = deepcopy(t1), deepcopy(t2)
    @inbounds for j in eachindex(dt), i in eachindex(phi)
        lam1[i,j], lam2[i,j] = compute_eigvals(t1, t2, phi[i], dt[j], T1, T2)
    end
    # Best parameters
    ip, idt = minloc2(lam2)
    # Find source polarisation and error therein
    spol, dspol = sourcepol(t1, t2, phi[ip], dt[idt])
    # Return minimum values as well
    (phi=phi, dt=dt, lam1=lam1, lam2=lam2, phi_best=phi[ip], dt_best=dt[idt], spol=spol, dspol=dspol)
end

function apply_split!(t1::SACtr, t2::SACtr, phi::Number, dt::Number)
    # Apply  a splitting operator to two SAC traces.  They are
    # assumed to be orthogonal, and the first trace should be anticlockwise
    # of the second.
    theta = phi - t1.cmpaz
    rotate_traces!(t1, t2, theta)
    # Move the fast trace forward
    shift_trace!(t1, -dt)
    rotate_traces!(t1, t2, -theta)
    nothing
end

function rotate_traces!(s1, s2, phi)
    # Version of SAC.rotate_through! which does not update headers
    phir = (eltype(s1.t))(deg2rad(phi))
    sinp, cosp = sincos(phir)
    for i = 1:s1.npts
        @inbounds s1.t[i], s2.t[i] = cosp*s1.t[i] - sinp*s2.t[i], sinp*s1.t[i] + cosp*s2.t[i]
    end
end

function shift_trace!(t, dt)
    # Version of SAC.tshift! which does not update headers and zeros out
    # points shifted from outside the window
    T = eltype(t.t)
    n = round(Int, dt/t.delta)
    if n == 0
        return s
    end
    _arrayshift!(t.t, n)
end

function _arrayshift!(a, n)
    # Shift the elements of an array back by n items
    N = length(a)
    O = zero(eltype(a))
    if n > 0
        @inbounds for i in N:-1:(n + 1)
            a[i] = a[i-n]
        end
        @inbounds for i in 1:n
            a[i] = O
        end
    elseif n < 0
        @inbounds for i in 1:(N + n)
            a[i] = a[i-n]
        end
        @inbounds for i in (N + n + 1):N
            a[i] = O
        end
    end
    a
end

function compute_eigvals(t1::SACtr, t2::SACtr, phi::Number, dt::Number, T1::SACtr, T2::SACtr)
    # Compute the two eigenvalues for traces for a specific fast orientation, phi,
    # and delay time, dt.
    # Note that we swap the delay time, because we are applying the inverse
    # splitting operator.
    copy_trace!(t1, T1)
    copy_trace!(t2, T2)
    apply_split!(T1, T2, phi, -dt)
    c = covar(T1, T2)
    eval = eigvals(c)
    maximum(eval), minimum(eval)
end

function sourcepol(t1::SACtr, t2::SACtr, phi::Number, dt::Number)
    # Compute the source polarisation and error for the best-fitting phi and dt
    T1, T2 = deepcopy(t1), deepcopy(t2)
    apply_split!(T1, T2, phi, -dt)
    c = covar(T1, T2)
    eval, evec = eigen(c)
    i1 = argmax(eval)
    i2 = 3 - i1
    lam1 = eval[i1]
    lam2 = eval[i2]
    spol = rad2deg(atan(evec[2,i1], evec[1,i1]))
    dspol = rad2deg(atan(lam2/lam1))
    spol, dspol
end

function covar(t1::SACtr, t2::SACtr)
    # Return the covariance matrix for two traces
    c11 = c22 = c12 = zero(eltype(t1.t))
    @inbounds for i in 1:t1.npts
        c11 += t1.t[i]^2
        c22 += t2.t[i]^2
        c12 += t1.t[i]*t2.t[i]
    end
    c = SArray{Tuple{2,2},SAC.SACFloat}(c11, c12, c12, c22)
    c
end

@eval begin
    """
        copy_trace!(a::SACtr, b::SACtr)

    Copy the SACtr trace `a` into `b` without allocating.
    """
    function copy_trace!(a::SACtr, b::SACtr)
        $([:(b.$h = a.$h) for h in fieldnames(SACtr)[1:end-1]]...)
        b.t .= a.t
        nothing
    end
end

function minloc2(A)
    # Return an 2-tuple which contains the location of the minimum
    # value of a 2-dimensional array.
    # Why isn't this implemented (in general) for n-arrays in Julia already?
    imin = jmin = 0
    min = zero(eltype(A))
    for j in 1:size(A, 2)
        for i in 1:size(A, 1)
            if i == j == 1
                imin = jmin = 1
                min = A[i,j]
            else
                if A[i,j] < min
                    imin = i
                    jmin = j
                    min = A[i,j]
                end
            end
        end
    end
    imin, jmin
end

end # module
