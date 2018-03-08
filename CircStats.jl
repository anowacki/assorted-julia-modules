__precompile__()

"""
module CircStats contains routines for investigating the statistic of circular
data.

By default, all input and output from routines is in radians, but in general passing
`true` as the last argument to a routine will change this to degrees.
"""
module CircStats

import Distributions

export
    cdist,
    cmean,
    cmedian,
    cresultant,
    cstd,
    cvariance,
    fit_vonMises,
    von_mises_cdf,
    von_mises_pdf,
    watson_U2

"""`cdist(a, b, degrees::Bool=false) -> angle`

Return the angular distance from `a` to `b` (b - a) in the forward
direction; hence if `b` is 'behind' `a`, `distance` will be negative.

Angles are confined to be the smaller possible angle, so are in the range
[-π:π], or [-180°:180°].

Angles are in radians, unless `degrees` == true.
"""
cdist(a, b) = (b%2pi - a%2pi + pi)%2pi - pi
cdist(a, b, degrees::Bool) = degrees ? (b%360. - a%360. + 180.)%360. - 180. : cdist(a, b)

"`cmean(a::Array, degrees::Bool=false) -> mean`

Return the mean angle from the set of angles in `a`.

Angles are in radians, unless `degrees` == true.
"
cmean(a) = atan2(sum(sin.(a)), sum(cos.(a)))
cmean(a, degrees::Bool) = degrees ? rad2deg(cmean(deg2rad.(a))) : cmean(a)

"`cmedian(a::Array, degrees::Bool=false, axial=false) -> median`

Return the median angle from the set of angles `a`.

Angles are in radians, unless `degrees` == true.
"
function cmedian(a, degrees::Bool=false, axial::Bool=false)
    medians = Vector{eltype(a)}()
    A = sort(a)
    n = length(A)
    degrees && (A .= deg2rad.(A))
    axial && (A .= 2A)
    # Compute trial bisectors, which are either the points themselves or adjacent means
    p = Array{Float64}(n)
    if iseven(n)
        for i = 1:n
            j = (i + 1 <= n) ? i + 1 : 1
            p[i] = cmean([A[i], A[j]])
        end
    else
        p[:] .= A[:]
    end
    # Try all possible diameters
    for i = 1:n
        # Count points on either side of diameter
        n_plus = sum(cdist.(p[i], A) .> 0)
        n_minus = sum(cdist.(p[i], A) .< 0)
        if n_plus == n_minus
            # Determine which side of circle is correct direction by counting the number
            # within pi/2 of each of the two opposing possible medians
            if sum(abs.(cdist.(p[i], A)) .<= pi/2) > sum(abs.(cdist.(p[i], A)) .> pi/2)
                push!(medians, A[i])
            else
                push!(medians, (A[i] + 2pi)%2pi - pi)
            end
        end
    end
    # If there is more than one median, take the mean thereof (Otieno & Anderson-Cook
    # (2003), Journal of Modern Applied Statistical Methods, 2(1), 168-176)
    median = if length(medians) > 1
        cmean(medians)
    elseif length(medians) == 1
        medians[1]
    else
        error("Zero medians found.  Are data axial but axial!=true?")
    end
    degrees && (median = rad2deg(median))
    axial && (median /= 2)
    median
end

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
    sqrt(sum(sin.(deg2rad.(a)))^2 + sum(cos.(deg2rad.(a)))^2)/length(a) :
    sqrt(sum(sin.(a))^2 + sum(cos.(a))^2)/length(a)
cresultant(a, w::Real, degrees::Bool=false) = degrees ?
    cresultant(a, true)*deg2rad(w)/(2sin(deg2rad(w)/2)) : cresultant(a)*w/(2sin(w/2))

"""
    cstd(a, degrees=false) -> σ

Return the standard deviation, `σ`, for a set of angles `a`.

Angles are in radians, unless `degrees == true`.
"""
cstd(a, degrees::Bool=false) = sqrt(-2.*log(cresultant(a, degrees)))

"""
    cvariance(a, degrees=false) -> σ²

Return the circular variance, `σ²`, of a set of angles `a`.

Angles are in radians, unless `degrees == true`.
"""
function cvariance(a, degrees::Bool=false)
    a_mean = cmean(a, degrees)
    if degrees
        1. - sum(cos.(deg2rad.(a - a_mean)))/length(a)
    else
        1. - sum(cos.(a - a_mean))/length(a)
    end
end

"
`von_mises_pdf(a, µ, κ, degrees=false) -> p`

Return the Von Mises probability density function at `a`, for a Von Mises distribution
with circular mean `µ` and concentration `κ`.

Angles are in radians, unless `degrees` == true.
"
von_mises_pdf(a, mu::Real, k::Real, degrees::Bool=false) = degrees ?
                                    exp(k*cos(deg2rad(a - mu)))/(2pi*besseli(0, k)) :
                                    exp(k*cos(a - mu))/(2pi*besseli(0, k))

"""
    von_mises_cdf(a, μ, κ, degrees=false) -> cdf

Return the von Mises cumulative distribution function at
`a` for a distribution with mean `μ` and concentration `κ`.
"""
function von_mises_cdf(a, μ, κ, degrees=false)
    degrees && (μ = deg2rad(μ))
    d = Distributions.VonMises(μ, κ)
    Distributions.cdf(d, degrees ? deg2rad(a) : a)
end

"""
    watson_U2(θ, cdf, α=0.05, degrees=false; axial=false) -> test, U², U²crit

Perform Watson's U²n test to determine if a sample of angles `θ`
fit the supplied cdf of a distribution at the `α` level of
significance.  `test` is `true` if the fit *is* significant.
Return also the value of the test statistic `U²` and the critical value
`U²crit`.

    watson_U2(θ, α=0.05, degrees=false; axial=false) -> test, U², U²crit

Perform the test against the maximum-likelihood von Mises distribution
for the sample data.  This is fit using `fit_vonMises()`.
"""
function watson_U2(θ, cdf::Function, α=0.05, degrees=false; axial=false)
    axial && (θ = 2θ)
    degrees && (θ = deg2rad.(θ))
    θ = sort(mod.(θ, 2π))
    n = length(θ)
    V = cdf.(θ) - cdf(0)
    V̄ = sum(V)/n
    U² = sum(V.^2) - sum((2*(1:n) - 1).*V/n) + n*(1/3 - (V̄ - 1/2)^2)
    U²crit = watson_U2_crit(n, α)
    U² > U²crit, U², U²crit
end

function watson_U2(θ, α::Number=0.05, degrees=false; kwargs...)
    μ, κ = fit_vonMises(θ, degrees; kwargs...)
    watson_U2(θ, x->von_mises_cdf(x, μ, κ, degrees), α, degrees; kwargs...)
end

"Table of critical values of Watson's U² test"
const WATSON_U2_TABLE = [0.143 0.000 0.161 0.164 0.165
                         0.145 0.173 0.194 0.213 0.224
                         0.146 0.176 0.202 0.233 0.252
                         0.148 0.177 0.205 0.238 0.262
                         0.149 0.179 0.208 0.243 0.269
                         0.149 0.180 0.210 0.247 0.274
                         0.150 0.181 0.211 0.250 0.278
                         0.150 0.182 0.212 0.252 0.281
                         0.150 0.182 0.213 0.254 0.283
                         0.150 0.183 0.215 0.256 0.287
                         0.151 0.184 0.216 0.258 0.290
                         0.151 0.184 0.216 0.259 0.291
                         0.151 0.184 0.217 0.259 0.292
                         0.151 0.185 0.217 0.261 0.293
                         0.152 0.185 0.219 0.263 0.296
                         0.152 0.186 0.219 0.264 0.298
                         0.152 0.186 0.220 0.265 0.299
                         0.152 0.186 0.221 0.266 0.301
                         0.152 0.187 0.221 0.267 0.302]
"Values of α (columns) for WATSON_U2_TABLE"
const WATSON_U2_ALPHAS = (0.100, 0.050, 0.025, 0.010, 0.005)
"Values of n (rows) for WATSON_U2_table"
const WATSON_U2_N = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 30, 40, 50, 100, 200)

"""
    watson_U2_crit(n, α) -> U²crit

Return the critical value at the α level for a sample of `n` angles.
"""
function watson_U2_crit(n, α)
    2 <= n || throw(ArgumentError("`n` (supplied $n) must be 2 or more"))
    any(isapprox.(α, WATSON_U2_ALPHAS, atol=0.0001)) == 1 ||
        throw(ArgumentError("Significance level for Watson's U² test must " *
                            "be one of $(WATSON_U2_ALPHAS); asked for $α"))
    indn = indmin(abs.(n .- WATSON_U2_N))
    indα = indmin(abs.(α .- WATSON_U2_ALPHAS))
    WATSON_U2_TABLE[indn,indα]
end

"""
    test_two_populations_differ(θ₁, θ₂, α, degrees=false; axial=false) -> ::Bool

Perform the Watson-Williams test to determine if two independent
circular measurements differ significantly from each other.

Assumes that samples are drawn from a von Mises distribution with a
concentration parameter κ that is the same, and κ > 2.
"""
function test_two_populations_differ(a, b, α, degrees=false; axial=false)
end

"""
    fit_vonMises(θ, degrees=false; axial=false) -> μ, κ

Return the maximum-likelihood values of the mean `μ` and concentration
`κ` for the von Mises distribution which fits the set of angles `\theta`.
"""
function fit_vonMises(θ, degrees=false; axial=false)
    degrees && (θ = deg2rad.(θ))
    axial && (θ = 2θ)
    μ = cmean(θ)
    κ = estimate_kappa_vonMises(θ, μ)
    degrees ? (rad2deg(μ), κ) : (μ, κ)
end

"""
    estimate_kappa_vonMises(θ, μ) -> κ

Return the maximum likelihood estimation of the concentration
parameter `κ` for the best-fitting von Mises distribution, given
a set of angles `θ` in radians.

Uses the polynomial approximation given by Best and Fisher (1981).
**N.B.** This may not be reliable when R̄ is small (e.g., < 0.7).
"""
function estimate_kappa_vonMises(θ, μ)
    R̄ = mean(cos.(θ - μ))
    if 0 <= R̄ < 0.53
        2R̄ + R̄^3 + 5R̄^5/6
    elseif 0.53 <= R̄ < 0.85
        -0.4 + 1.39R̄ + 0.43/(1 - R̄)
    elseif R̄ >= 0.85
        1/(R̄^3 - 4R̄^2 + 3R̄)
    end
end

end # module
