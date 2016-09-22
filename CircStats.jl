"""
module CircStats contains routines for investigating the statistic of circular
data.

By default, all input and output from routines is in radians, but in general passing
`true` as the last argument to a routine will change this to degrees.
"""
module CircStats

__precompile__()

export
    cdist,
    cmean,
    cmedian,
    cresultant,
    von_mises_pdf

"""`cdist(a, b, degrees::Bool=false) -> angle`

Return the angular distance from `a` to `b` (b - a) in the forward
direction; hence if `b` is 'behind' `a`, `distance` will be negative.

Angles are confined to be the smaller possible angle, so are in the range
[-π:π], or [-180°:180°].

Angles are in radians, unless `degrees` == true.
"""
cdist(a, b) = (b%2pi - a%2pi + pi)%2pi - pi
cdist(a, b, degrees::Bool) = degrees ? (b%360 - a%360 + 180)%360 - 180 : cdist(a, b)

"`cmean(a::Array, degrees::Bool=false) -> mean`

Return the mean angle from the set of angles in `a`.

Angles are in radians, unless `degrees` == true.
"
cmean(a) = atan2(sum(sin(a)), sum(cos(a)))
cmean(a, degrees::Bool) = degrees ? rad2deg(cmean(deg2rad(a))) : cmean(a)

"`cmedian(a::Array, degrees::Bool=false) -> median`

Return the median angle from the set of angles `a`.

Angles are in radians, unless `degrees` == true.
"
function cmedian(a)
    medians = Float64[]
    A = sort(a)
    n = length(A)
    # Compute trial bisectors, which are either the points themselves or adjacent means
    p = Array{Float64}(n)
    if iseven(n)
        for i = 1:n
            j = (i + 1 <= n) ? i + 1 : 1
            p[i] = cmean([A[i], A[j]])
        end
    else
        p[:] = A[:]
    end
    # Try all possible diameters
    for i = 1:n
        # Count points on either side of diameter
        n_plus = sum(cdist(p[i], A) .> 0)
        n_minus = sum(cdist(p[i], A) .< 0)
        if n_plus == n_minus
            # Determine which side of circle is correct direction by counting the number
            # within pi/2 of each of the two opposing possible medians
            if sum(abs(cdist(p[i], A)) .<= pi/2) > sum(abs(cdist(p[i], A)) .> pi/2)
                push!(medians, A[i])
            else
                push!(medians, (A[i] + 2pi)%2pi - pi)
            end
        end
    end
    # If there is more than one median, take the mean thereof (Otieno & Anderson-Cook
    # (2003), Journal of Modern Applied Statistical Methods, 2(1), 168-176)
    if length(medians) > 1
        cmean(medians)
    else
        medians[1]
    end
end
cmedian(a, degrees::Bool) = degrees ? rad2deg(cmedian(deg2rad(a))) : cmedian(a)

"Example sorted angle data, quoted in Otieno & Anderson-Cook (2003)"
frog_data = Float64[104, 110, 117, 121, 127, 130, 136, 144, 152, 178, 184, 192, 200, 316]


function test_median(;plot::Bool=false)
    # Test for median computation from PyCircStat:
    # https://github.com/circstat/pycircstat/blob/master/tests/test_descriptive.py
    alpha = [3.73153000, 1.63904879, 4.03175622, 3.90422402, 4.61029613,
             4.04117818, 5.79313473, 5.50863002, 5.81530225, 2.44973903,
             2.12868554, 0.09073566, 0.05581025, 5.10673712, 1.68712454,
             3.72915575, 4.45439608, 4.70694685, 3.58470730, 2.49742028]
    mall = mod1(-2.2467, 2pi)
    median = cmedian(alpha)[1]
    if plot
        cplot_pts(alpha)
        PyPlot.plot(cos(median), sin(median), "o", markersize=20, label="CircStats")
        PyPlot.plot(cos(mall), sin(mall), "o", markersize=20, label="PyCircStat test")
        PyPlot.legend()
    end
    @assert isapprox(median, mall, atol=1e-3) "$(mall) != $(median)"

end

"
`cresultant(a, degrees=false) -> R`
`cresultant(a, w, degrees=false) -> Rc`

Return the resultant vector length, `R`, from a set of angles, `a`.

If data are binned by binwidth `w`, an unbiased estimate of `R`, `Rc` is returned
when `w` is supplied.

Angles are in radians, unless `degrees` == true.
"
cresultant(a, degrees::Bool=false) = degrees ?
    sqrt(sum(sin(deg2rad(a)))^2 + sum(cos(deg2rad(a)))^2)/length(a) :
    sqrt(sum(sin(a))^2 + sum(cos(a))^2)/length(a)
cresultant(a, w::Number, degrees::Bool=false) = degrees ?
    cresultant(a, true)*deg2rad(w)/(2sin(deg2rad(w)/2)) : cresultant(a)*w/(2sin(w/2))

"
`von_mises_pdf(a, µ, κ, degrees=false) -> p`

Return the Von Mises probability density function at `a`, for a Von Mises distribution
with circular mean `µ` and concentration `κ`.

Angles are in radians, unless `degrees` == true.
"
von_mises_pdf(a, mu::Number, k::Number, degrees::Bool=true) = degrees ?
                                    exp(k*cos(deg2rad(a - mu)))/(2pi*besseli(0, k)) :
                                    exp(k*cos(a - mu))/(2pi*besseli(0, k))


end # module
