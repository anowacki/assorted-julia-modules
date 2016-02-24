"""
module CircPlot provides a number of ways of plotting circular data.

By default, all input to routines is in radians, but in general passing
`true` as the last argument to a routine will change this to degrees.
"""
module CircPlot

import PyPlot

export
    cplot_histogram,
    cplot_pts

"
`cplot_pts(a::Array, degrees=false; azimuth=false, hold=false, args...)`

Use PyPlot to plot a set of points on the circumference of the unit circle.

Use `azimuth=true` to plot angles clockwise from north.

Use `hold=true` to overlay additional plots on an existing plot.

Supply any additional arguments to `PyPlot.plot()` as keyword arguments; for example,
to set the size of points, do

    cplot_pts(angles; markersize=20)
"
function cplot_pts(a, degrees::Bool=false;
                   azimuth::Bool=false, hold::Bool=false, args...)
    hold || PyPlot.clf()
    if degrees
        xfunc = cosd
        yfunc = sind
        circum = 0:1:361
    else
        xfunc = cos
        yfunc = sin
        circum = 0:pi/100:2pi+pi/100
    end
    if azimuth xfunc, yfunc = yfunc, xfunc end
    hold || PyPlot.plot(xfunc(circum), yfunc(circum))
    PyPlot.plot(xfunc(a), yfunc(a), "o"; args...)
    if !hold
        PyPlot.axis("square")
        PyPlot.xlim(-1.05, 1.05)
        PyPlot.ylim(-1.05, 1.05)
        if azimuth
            PyPlot.xlabel("East")
            PyPlot.ylabel("North")
        else
            PyPlot.xlabel("x")
            PyPlot.ylabel("y")
        end
    end
    return
end

"
`cplot_histogram(a, w, degrees=false; azimuth=false, axial=false, args...)`

Plot a polar histogram of a set of angles `a` using a binwidth `w`.  If 2π (or 360°)
is not an integer multiple of `w`, then the nearest value for which this is true is
selected (and a warning issued).

Use `azimuth=true` to plot angles clockwise from north.

Use `axial=true` if data are π-ambiguous (i.e., are orientational rather than
    directional data).

Supply any additional arguments to `PyPlot.plot()` as keyword arguments.
"
function cplot_histogram(a, w, degrees::Bool=false;
                         azimuth::Bool=false, axial::Bool=false, bottom=0, args...)
    n = (degrees ? round(Int, 360/w) : round(Int, 2pi/w)) + 1
    w = 2pi/(n - 1)
    bins = linspace(0, 2pi, n)
    data = degrees ? deg2rad(mod(a, 360)) : mod(a, 2pi)
    data = axial ? [data; mod(data + pi, 2pi)] : data
    bins, count = hist(data, bins)
    PyPlot.clf()
    ax = PyPlot.plt[:subplot](111, polar=true)
    if azimuth
        ax[:set_theta_zero_location]("N")
        ax[:set_theta_direction](-1)
    end
    bars = ax[:bar](bins[1:end-1], count, width=w, bottom=bottom)
    PyPlot.plt[:show]()
    return
end

end # module
