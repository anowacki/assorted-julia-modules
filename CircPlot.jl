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

const piby2 = π/2

#=
Exported routines
=#
"""
    cplot_pts(a::Array, degrees=false; azimuth=false, hold=false, args...) -> ::Plots.Plot

Use Plots to plot a set of points on the circumference of the unit circle.

Use `azimuth=true` to plot angles clockwise from north.

Use `hold=true` to overlay additional plots on an existing plot.

Supply any additional arguments to `PyPlot.plot()` as keyword arguments; for example,
to set the size of points, do

    julia> cplot_pts(angles; markersize=20)
"""
function cplot_pts(a, degrees::Bool=false;
                   azimuth::Bool=false, axial::Bool=false, line::Bool=true,
                   padding=0.05, R=1, centre=(0,0), kwargs...)
    circum = 0:0.01:2π
    markerstyle = (0.1,:blue)
    degrees && (a = deg2rad(a))
    azimuth && (a = piby2 - a)
    R += padding
    x0, y0 = centre
    p = plot(aspect_ratio=:equal,
        xlim=(-R+x0,R+x0), ylim=(-R+y0,R+y0), legend=false)
    line && plot!(p, x0+cos(circum), y0+sin(circum), l=:black)
    scatter!(p, x0+cos(a), y0+sin(a), m=markerstyle; kwargs...)
    axial && scatter!(p, x0+cos(a+π), y0+sin(a+π), m=markerstyle; kwargs...)
    if azimuth
        xlabel!(p, "East")
        ylabel!(p, "North")
    else
        xlabel!(p, "x")
        ylabel!(p, "y")
    end
    p
end

"
    cplot_histogram(a, binwidth, degrees=false; azimuth=false, axial=false; kwargs...) -> ::Plots.Plot

Plot a polar histogram of a set of angles `a` using bins `binwidth` wide.
`binwidth` should be given in the same convention (radians or °) as the data.
If π (or 180°) is not an integer multiple of `binwidth`, then the nearest value
for which this is true is selected.

Use `azimuth=true` to plot angles clockwise from north.

Use `axial=true` if data are π-ambiguous (i.e., are orientational rather than
directional data).

Set options for the surrounding circle using `circ=Dict(:opt=>val)`, where
the keys and values in the Dict will be passed to Plots for plotting the circle
shape around the plot.

Additional keyword arguments `kwargs...` are passed to `Plots.plot()` when
plotting the bars.
"
function cplot_histogram(a, binwidth, degrees::Bool=false;
                         azimuth::Bool=false, axial::Bool=false,
                         circ=nothing,
                         kwargs...)
    hist = fit_hist(a, binwidth, degrees, azimuth, axial)
    n = length(hist.weights)
    edges = hist.edges[1]
    R = maximum(hist.weights)
    p = plot(aspect_ratio=:equal, xlim=(-R,R), ylim=(-R,R), legend=false)
    circ != nothing && plot!(p, circle(R); circ...)
    plot!(p, sector.(edges[1:n], edges[2:n+1], hist.weights), l=stroke(0),
        fill=:black; kwargs...)
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
the radius (in plot units) of the maximum bin.  Set `minobs` to be the minimum
bar count required for the plot to fill the entire maximum radius from `maxr`.
"""
function chistogram(a, binwidth, degrees=false, x=0, y=0;
        axial=false, azimuth=false, scale=1, maxr=nothing, minobs=1)
    hist = fit_hist(a, binwidth, degrees, azimuth, axial)
    n = length(hist.weights)
    edges = hist.edges[1]
    R = maximum(hist.weights)
    maxr != nothing && (scale = maxr/max(minobs, R))
    sector.(edges[1:n], edges[2:n+1], hist.weights*scale, x, y)
end

#=
Unexported routines
=#
"""Return a `StatsBase.Histogram` with the results of binning the angles `a`
into bins `binwidth` wide."""
function fit_hist(a, binwidth, degrees::Bool, azimuth::Bool, axial::Bool)
    n = (degrees ? round(Int, 360/binwidth) : round(Int, 2pi/w))
    binwidth = 2pi/n
    bins = linspace(0, 2pi, n+1)
    data = degrees ? deg2rad.(mod.(a, 360)) : mod.(a, 2pi)
    axial && (data = [mod.(data, pi); mod.(data, pi).+pi])
    azimuth && (data .= mod.(piby2 .- data, 2pi))
    StatsBase.fit(StatsBase.Histogram, data, bins, closed=:left)
end

"""Return a sensible number of segments between two angles in radians θ₁ and θ₂."""
number_of_segments(θ₁, θ₂) = max(2, round(Int, abs(θ₂ - θ₁)*90/π))

"""Return a circle shape with radius `r` centred at `x` and `y`."""
circle(r, x=0, y=0) = Shape(arc(0, 2π, r, x, y))

"""Return an arc (part of a circle's circumference) between `θ₁` and `θ₂` radians
anticlockwise from the +x axis, with radius `r`, centred at `x` and `y`."""
arc(θ₁, θ₂, r=1, x=0, y=0) = Tuple{Float64,Float64}[(x+r*cos(θ), y+r*sin(θ))
    for θ in linspace(θ₁, θ₂, number_of_segments(θ₁, θ₂))]

"""Return a sector (pie wedge) between angles `θ₁` and `θ₂` radians from x with
radius `r`, centred at `x` and `y`."""
sector(θ₁, θ₂, r=1, x=0, y=0) = Shape(vcat(arc(θ₁, θ₂, r, x, y), (x, y)))

"""Return a filled sector (chordshape) between two angles θ₁ and θ₂ radians from
x and two radii, `r₁` and `r₂`, centred at `x` and `y`."""
chordshape(θ₁, θ₂, r₁, r₂, x=0, y=0) =
    Shape(vcat(arc(θ₁, θ₂, r₁, x, y), arc(θ₂, θ₁, r₂, x, y)))

end # module
