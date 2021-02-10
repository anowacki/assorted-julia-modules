"""
## CircPlot

Plot circular data in a number of ways:

- `cplot_histogram`: Plot a polar histogram of angular data.
- `cplot_pts`: Plot angular data as points on the edge of a circle.
- `chistogram`: Get a set of `Plots.Shape`s which can be plotted as a polar
  histogram at any point on a map, for example.

By default, all input to routines is in radians, but all routines have a positional
argument to set `degrees` to `true` when not using radians.
"""
module CircPlot

import StatsBase

using Plots

export
    chistogram,
    cplot_histogram,
    cplot_pts

#=
Exported routines
=#
"""
    cplot_pts(a::Array, degrees=false; azimuth=false, hold=false, args...) -> ::Plots.Plot

Use Plots to plot a set of points on the circumference of the unit circle.

Use `azimuth=true` to plot angles clockwise from north.

Supply any additional arguments to `Plots.plot()` as keyword arguments; for example,
to set the size of points, do

    julia> cplot_pts(angles; markersize=20)
"""
function cplot_pts(a, degrees::Bool=false;
                   azimuth::Bool=false, axial::Bool=false, line::Bool=true,
                   padding=0.05, R=1, centre=(0,0), kwargs...)
    circum = 0:0.01:2π
    markerstyle = (0.1,:blue)
    degrees && (a = deg2rad.(a))
    azimuth && (a = π/2 .- a)
    R += padding
    x0, y0 = centre
    p = plot(aspect_ratio=:equal,
        xlim=(-R+x0,R+x0), ylim=(-R+y0,R+y0), legend=false)
    line && plot!(p, x0.+cos.(circum), y0.+sin.(circum), l=:black)
    scatter!(p, x0.+cos.(a), y0.+sin.(a), m=markerstyle; kwargs...)
    axial && scatter!(p, x0.+cos.(a+π), y0.+sin.(a+π), m=markerstyle; kwargs...)
    if azimuth
        xlabel!(p, "East")
        ylabel!(p, "North")
    else
        xlabel!(p, "x")
        ylabel!(p, "y")
    end
    p
end

"""
    cplot_histogram(a, binwidth, degrees=false; azimuth=false, axial=false, weights=ones(length(a)), circ=nothing; kwargs...) -> ::Plots.Plot

Plot a polar histogram of a set of angles `a` using bins `binwidth` wide.
`binwidth` should be given in the same convention (radians or °) as the data.
If π (or 180°) is not an integer multiple of `binwidth`, then the nearest value
for which this is true is selected.

Use `azimuth=true` to plot angles clockwise from north.

Use `axial=true` if data are π-ambiguous (i.e., are orientational rather than
directional data).

`weights` must have the same length as `a`, and if passed in will scale each
point `i` by `weights[i]`.

Set options for the surrounding circle using `circ=Dict(:opt=>val)`, where
the keys and values in the Dict will be passed to Plots for plotting the circle
shape around the plot.

Additional keyword arguments `kwargs...` are passed to `Plots.plot()` when
plotting the bars.
"""
function cplot_histogram(a, binwidth, degrees::Bool=false;
                         azimuth::Bool=false, axial::Bool=false,
                         weights=ones(eltype(a), length(a)), circ=nothing,
                         kwargs...)
    hist = fit_hist(a, binwidth, degrees, axial, weights=weights)
    n = length(hist.weights)
    edges = hist.edges[1]
    azimuth && (edges = π/2 .- collect(edges))
    R = maximum(hist.weights)
    p = plot(aspect_ratio=:equal, xlim=(-R,R), ylim=(-R,R), legend=false)
    circ != nothing && plot!(p, circle(R); circ...)
    plot!(p, sector.(edges[1:n], edges[2:n+1], hist.weights), fill=:black; kwargs...)
    p
end

"""
    chistogram(a, binwidth, degrees::Bool=true, x=0, y=0; axial=false, azimuth=false, scale=1, maxr=nothing) -> bars::Array{Shape}

Return a set of Shapes which show a polar histogram created from the set of
angles `a` with bins `binwidth` wide, centred at `x` and `y`.  These shapes can
be plotted with `plot(bars)`.

Set `axial` to `true` for data with π-ambiguity and `azimuth` to `true` if data are
clockwise from north.

`scale` determine what length to scale the radii by, or alternatively set `maxr` to
the radius (in plot units) of the maximum bin.  Set `minweight` to be the minimum
bar count required for the plot to fill the entire maximum radius from `maxr`.
"""
function chistogram(a, binwidth, degrees=false, x=0, y=0;
        weights=ones(a), axial=false, azimuth=false, scale=1, maxr=nothing,
        minweight=1, minobs=nothing)
    if minobs != nothing
        warning("`minobs` is deprecated: use `minweight` instead")
        minweight = minobs
    end
    hist = fit_hist(a, binwidth, degrees, axial, weights=weights)
    n = length(hist.weights)
    edges = hist.edges[1]
    R = maximum(hist.weights)
    maxr != nothing && (scale = maxr/max(minweight, R))
    azimuth && (edges = π/2 .- collect(edges))
    sector.(edges[1:n], edges[2:n+1], hist.weights*scale, x, y)
end

#=
Unexported routines
=#
"""    fit_hist(a, binwidth, degrees=false, axial=false; weights=ones(a))

Return a `StatsBase.Histogram` with the results of binning the angles `a`
into bins `binwidth` wide.  If `axial` is `true`, then directions have
π symmetry.  Optionally sypply a `weights` array containing the weighting
for each value in `a`."""
function fit_hist(a, binwidth, degrees::Bool, axial::Bool;
                  weights=ones(eltype(a), length(a)))
    binwidth isa Bool &&
        throw(ArgumentError("`binwidth` supplied as a `Bool`, but should be a real number"))
    length(weights) == length(a) ||
        throw(ArgumentError("`weights` must have same length as number of angles"))
    n = (degrees ? round(Int, 360/binwidth) : round(Int, 2pi/binwidth))
    n = max(1, n)
    binwidth = 2pi/n
    bins = range(0, stop=2pi, length=n+1)
    data = degrees ? deg2rad.(mod.(a, 360)) : mod.(a, 2pi)
    # FIXME: Many angles are given as whole numbers: this cludge helps
    # an issue where floating point values sometimes are binned one way,
    # sometimes another, giving polar histograms which depend on the range
    # in which angles are given.
    data .+= 10*eps(float(eltype(a)))
    axial && (data .= mod.(data, pi))
    h = StatsBase.fit(StatsBase.Histogram, data,
                      StatsBase.FrequencyWeights(weights), bins, closed=:left)
    axial && (h.weights[end÷2+1:end] .= h.weights[1:end÷2])
    h
end

"""Return a sensible number of segments between two angles in radians θ₁ and θ₂."""
number_of_segments(θ₁, θ₂) = max(2, round(Int, abs(θ₂ - θ₁)*90/π))

"""Return a circle shape with radius `r` centred at `x` and `y`."""
circle(r, x=0, y=0) = Shape(arc(0, 2π, r, x, y))

"""Return an arc (part of a circle's circumference) between `θ₁` and `θ₂` radians
anticlockwise from the +x axis, with radius `r`, centred at `x` and `y`."""
arc(θ₁, θ₂, r=1, x=0, y=0) = Tuple{Float64,Float64}[(x+r*cos(θ), y+r*sin(θ))
    for θ in range(θ₁, stop=θ₂, length=number_of_segments(θ₁, θ₂))]

"""Return a sector (pie wedge) between angles `θ₁` and `θ₂` radians anticlockwise
from the x-axis with radius `r`, centred at `x` and `y`."""
sector(θ₁, θ₂, r=1, x=0, y=0)::Shape =
    Shape(vcat(arc(θ₁, θ₂, r, x, y), Tuple{Float64,Float64}((x, y))))

"""Return a filled sector (arcshape) between two angles θ₁ and θ₂ radians
anticlockwise from the x-axis, and two radii, `r₁` and `r₂`, centred at `x` and `y`."""
arcshape(θ₁, θ₂, r₁, r₂, x=0, y=0) =
    Shape(vcat(arc(θ₁, θ₂, r₁, x, y), arc(θ₂, θ₁, r₂, x, y)))

end # module
