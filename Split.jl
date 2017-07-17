module Split
# Perform shear wave splitting analysis on two orthogonal SAC traces

using SAC

# Module defaults
const split_ndt = 51
const split_nphi = 181
const split_dt_max = 4.

function search_phi(t1::SACtr, t2::SACtr;
                    nphi=split_nphi, ndt=split_ndt, dt_max=split_dt_max)
    # Return the array of minimum eigenvalues of the covariance matrices
    # when examining trial operators with nphi fast orientations
    phi = linspace(-90., 90., nphi)
    dt = zeros(ndt)
    lam1, lam2 = zeros(nphi,ndt), zeros(nphi,ndt)
    T1, T2 = deepcopy(t1), deepcopy(t2)
    for i = 1:nphi
        dt, lam1[i,:], lam2[i,:] = search_dt(t1, t2, phi[i], dt_max, ndt, T1, T2)
    end
    # Best parameters
    ip, idt = minloc2(lam2)
    # Find source polarisation and error therein
    spol, dspol = sourcepol(t1, t2, phi[ip], dt[idt])
    # Return minimum values as well
    return phi, dt, lam1, lam2, phi[ip], dt[idt], spol, dspol
end

function search_dt(t1::SACtr, t2::SACtr, phi, dt_max, ndt, T1::SACtr, T2::SACtr)
    # For a given fast orientation, return the smaller eigenvalue of the
    # covariance matrix in an array accompanied by a vector of dts.
    dt = linspace(dt_max/ndt, dt_max, ndt)
    lam1, lam2 = zeros(ndt), zeros(ndt)
    for i = 1:ndt
        lam1[i], lam2[i] = compute_eigvals(t1, t2, phi, dt[i], T1, T2)
    end
    return dt, lam1, lam2
end

function apply_split!(t1::SACtr, t2::SACtr, phi::Number, dt::Number)
    # Apply  a splitting operator to two SAC traces.  They are
    # assumed to be orthogonal, and the first trace should be anticlockwise
    # of the second.
    theta = phi - t1.cmpaz
    SAC.rotate_through!(t1, t2, theta)
    # Move the fast trace forward
    SAC.tshift!(t1, -dt)
    SAC.rotate_through!(t1, t2, -theta)
    return
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
    eval = Base.eigvals(c)
    return maximum(eval), minimum(eval)
end

function sourcepol(t1::SACtr, t2::SACtr, phi::Number, dt::Number)
    # Compute the source polarisation and error for the best-fitting phi and dt
    T1, T2 = deepcopy(t1), deepcopy(t2)
    apply_split!(T1, T2, phi, -dt)
    c = covar(T1, T2)
    eval, evec = Base.eig(c)
    i1 = indmax(eval)
    i2 = 3 - i1
    lam1 = eval[i1]
    lam2 = eval[i2]
    spol = rad2deg(atan2(evec[2,i1], evec[1,i1]))
    dspol = rad2deg(atan(lam2/lam1))
    spol, dspol
end

function covar(t1::SACtr, t2::SACtr)
    # Return the covariance matrix for two traces
    c = Array{SAC.SACFloat}(2, 2)
    c[1,1] = sum(t1.t.^2)
    c[2,2] = sum(t2.t.^2)
    c[1,2] = c[2,1] = sum(t1.t.*t2.t)
    return c
end

const headers = "a"
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
    for j = 1:size(A, 2)
        for i = 1:size(A, 1)
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
    return (imin,jmin)
end

end # module
