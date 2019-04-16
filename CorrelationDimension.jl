module CorrelationDimension

import Plots

using Statistics: mean, covm, varm

using Distributions: fit
using LinearAlgebra: norm
using StatsBase: Histogram

export correlation_dimension, plot_correlation

"""
    correlation_dimension(x; dist=(x1,x2)->norm(x1.-x2)) -> (corrdim, logr, logNr, intercept)

Compute an estimate of the correlation dimension for a set of points in space, with
coordinates for the ith point given as `x[:,i]`.

This is given by

δ̂ = ∂[log N(r)]/∂[log r]

where N(r) is the number of point pairs within a distance r of each other.  This is
estimated by binning the pair-distances into bins (chosen automatically)
and fitting a straight line through the log-distribution.

The distance function `dist` by default is the Cartesian 2-norm, but any function may be
supplied.

The function returns a named tuple of the estimate of the correlation dimension, `corrdim`,
the lower edges of the distance bins, `logr`, the weights in the distance bins `logNr`, and
the `intercept` of the best fitting line (of which the slope is `corrdim`).
"""
function correlation_dimension(x; dist=(a,b)->norm(a .- b))
    D, n = size(x)
    # List of log(pairwise distances)
    point1 = Vector{float(eltype(x))}(undef, D)
    point2 = similar(point1)
    logr = float(eltype(x))[]
    for i in 1:n, j in (i+1):n
        point1 .= x[:,i]
        point2 .= x[:,j]
        push!(logr, log(dist(point1, point2)))
    end
    # Construct histogram
    h = fit(Histogram, logr)
    # Linear regression
    logr = h.edges[1][1:end-1]
    logNr = log.(cumsum(h.weights./(n*(n-1)/2)))
    c, δ̂ = linear_regression(logr, logNr)
    (corrdim=δ̂, logr=logr, logNr=logNr, intercept=c)
end

"""
    correlation_dimension(x1, x2, ..., xN; kwargs...) -> δ̂

Specify coordinates with dimension N as N separate vectors.
"""
correlation_dimension(x, y, args...; kwargs...) = correlation_dimension(permutedims(hcat(x, y, args...)))

"""
    plot_correlation((logr, logNr, slope, intercept))

Create a plot of the output from `correlation_dimension`.  Input is a tuple
of the log of the distance bins (`logr`), the weights of the histogram (`logNr`),
plus the `slope` and `intercept` of the best fitting line.  The `slope` therefore
is the estimate of the ccorrelation dimension.
"""
function plot_correlation((corrdim, logr, logNr, intercept))
    Plots.scatter(logr, logNr, label="Data", legend=:topleft)
    Plots.plot!(x->corrdim*x+intercept, label="y = $(round(corrdim, sigdigits=3))x + $(round(intercept, sigdigits=3))")
end

"""
    linear_regression(x, y)

Perform simple linear regression using Ordinary Least Squares. Returns `a` and `b` such
that `a + b*x` is the closest straight line to the given points `(x, y)`, i.e., such that
the squared error between `y` and `a + b*x` is minimized.
"""
function linear_regression(x::AbstractVector, y::AbstractVector)
    size(x) == size(y) || throw(DimensionMismatch("x and y must be the same size"))
    mx, my = mean(x), mean(y)
    b = covm(x, mx, y, my)/varm(x, mx)
    a = my - b*mx
    return a, b
end

end # module
