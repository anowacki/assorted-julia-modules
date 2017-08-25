__precompile__()

module SWSs

import Base: filter!, filter, getindex
using Base.Dates

using AutoHashEquals
using Geodesy

import SphericalGeom

export
    AbstractSWS,
    SWS,
    SWSTime,
    azimuth,
    backazimuth,
    midpoint,
    read_sheba,
    read_sheba_time,
    surf_dist,
    write_sheba


#=
    Types
=#
abstract type AbstractSWS end

"""
    SWS

Type containing information about a shear wave splitting measurement.
"""
@auto_hash_equals immutable SWS <: AbstractSWS
    evlo::Float64 # Degrees
    evla::Float64
    evdp::Float64 # km below sea level
    stlo::Float64
    stla::Float64
    stdp::Float64 # km below sea level
    stnm::String
    phi::Float64 # degrees
    dphi::Float64
    dt::Float64 # s
    ddt::Float64
    spol::Float64
    dspol::Float64
    ex::Float64 # m relative to
    ey::Float64
    ez::Float64
    sx::Float64
    sy::Float64
    sz::Float64
    dist::Float64 # Distance in m from source to receiver, straight line
    incidence::Float64 # Degrees from the downward direction
end

"""
    SWS(evlo, evla, evdp, stlo, stla, stdp, stnm, phi, dphi, dt, ddt, spol, dspol, origin) -> ::SWSobs

Return an `SWS` instance with extra fields filled in.
`origin` is a Geodesy.LLA instance and must be given
"""
function SWS(evlo, evla, evdp, stlo, stla, stdp, stnm, phi, dphi, dt, ddt, spol, dspol, origin)
    e = ENU(LLA(evla, evlo, -evdp*1000.), origin, wgs84)
    s = ENU(LLA(stla, stlo, -stdp*1000.), origin, wgs84)
    dist = distance(e, s)
    incidence = incidence_angle(e, s)
    SWS(evlo, evla, evdp, stlo, stla, stdp, stnm, phi, dphi, dt, ddt, spol, dspol,
           e.e, e.n, e.u, s.e, s.n, s.u, dist, incidence)
end

"Allow indexing into arrays of SWS by using the symbol form of the fieldname"
getindex(A::Array{SWS}, s::Symbol) = Array{typeof(getfield(A[1], s))}([getfield(a, s) for a in A])
getindex(t::SWS, s::Symbol) = getfield(t, s)
"Allow indexing into arrays of SWS obs by using a boolean array"
function getindex(A::Array{SWS}, b::BitArray)
    @assert length(b) == length(A) "Indexing array does not have correct dimensions"
    out = Array{SWS,1}()
    for i in 1:length(A)
        if b[i] push!(out, A[i]) end
    end
    out
end

"""
    SWSTime

Type containing a shear wave splitting measurement as well as the date and time of
the event.

Access the usual SWS field with `.sws` and the date with `.time`.
"""
immutable SWSTime <: AbstractSWS
    sws::SWS
    time::DateTime
end

"""
    getindex(s::Union{SWSTime,Array{SWSTime}}, k::Symbol) -> v

Return an array of values `v` for the fields named `k` in an array of `SWSTime`s.

Use `:time` for `k` to get the `DateTime`s for each measurement; otherwise, use `:sws`
to get an array of `SWS`s; or finally use any fieldname in the `SWS` type to get that.
"""
getindex(s::Array{SWSTime}, k::Symbol) = k == :time ? [ss.time for ss in s] :
                                             k == :sws ? [ss.sws for ss in s] :
                                             [getfield(ss.sws, k) for ss in s]
getindex(s::SWSTime, k::Symbol) = k == :time ? s.time :
                                      k == :sws ? s.sws : getfield(s.sws, k)


#=
    Functions on AbstractSWSs
=#

"Return the azimuth from the station to the receiver"
azimuth(s::AbstractSWS) = SphericalGeom.azimuth(s[:evlo], s[:evla], s[:stlo], s[:stla], true)
"Return the backazimuth from the receiver to the station"
backazimuth(s::AbstractSWS) = SphericalGeom.azimuth(s[:stlo], s[:stla], s[:evlo], s[:evla], true)

"Return the surface distance between the event and station"
surf_dist(s::Union{AbstractSWS,Array{<:AbstractSWS}}) = sqrt.((s[:ex] .- s[:sx]).^2 + (s[:ey] .- s[:sy]).^2)
"Return the incidence angle, measured away from downwards, from the event to the station"
incidence_angle(e::ENU, s::ENU) = rad2deg(atan2(s.u - e.u, sqrt((s.e - e.e).^2 + (s.n - e.n).^2)))

"""
    filter!(a::Array{SWS}, args...; kwargs...) -> a
Remove from the array of SWS measurements `a` any measurements which do not match
the filters specified by the `args` and `kwargs`.

    filter(a::Array{SWS}, args...; kwargs...) -> a_copied

Return a copy of `a` without non-matching measurements.


### Arguments

Arguments must be tuples containing a function and two values, the upper and lower limits.
Function must accept one argument (a SWS) and return one value.

```julia
julia> horizontal_dist(s::SWS) = sqrt((s.ex - s.sx)^2 + (s.ey - s.sy)^2);

julia> filter!(splits, (horizontal_dist, 0, 20_000)); # Measurements with epicentral distance less than 20 km
```

### Keyword arguments

Each keyword argument should corresond to
a field in the SWS type, and to it should be passed a tuple of two values, the lower
and upper bound (inclusive) that field may take.

The exception is the String field `stnm`, which accepts only a regular expression.

```julia
julia> filter!(splits, evdp=(-Inf,100)); # Keep only events shallower than 100 km
```
"""
function filter!(a::Array{<:AbstractSWS}, args...; kwargs...)
    # args
    for arg in args
        length(arg) == 3 && typeof(arg[1]) <: Function ||
            throw(ArgumentError("Arguments must be length-3 tuples"))
    end
    # kwargs
    field_types = Dict{Symbol,DataType}()
    for (f,t) in zip(fieldnames(SWS), SWS.types)
        field_types[f] = t
    end
    for (k,v) in kwargs
        k in keys(field_types) && field_types[k] == String &&
            !(typeof(v) <: Union{Regex,String}) &&
            throw(ArgumentError("Value for $k must be a String or Regex"))
    end
    delete_indices = Int[]
    fields_done = false
    for (i,s) in enumerate(a)
        fields_done = false
        for (f,v1,v2) in args
            fields_done && continue
            if !(v1 <= f(s) <= v2)
                push!(delete_indices, i)
                fields_done = true
            end
        end
        for (k,v) in kwargs
            fields_done && continue
            if haskey(field_types, k) && field_types[k] == String
                if !ismatch(Regex(v), s[k])
                    push!(delete_indices, i)
                    fields_done = true
                end
            else
                if !(v[1] <= s[k] <= v[end])
                    push!(delete_indices, i)
                    fields_done = true
                end
            end
        end
    end
    deleteat!(a, delete_indices)
    a
end

filter(a::Array{<:AbstractSWS}, args...; kwargs...) = filter!(deepcopy(a), args...; kwargs...)

@doc (@doc filter!) filter

"""
    midpoint(s::Union{SWS,Array{SWS}}) -> x, y, z

Return the midpoint of a straight line between the source and receiver in cartesian coordinates.
"""
midpoint(s::Union{AbstractSWS,Array{<:AbstractSWS}}) =
    (s[:sx] + s[:ex])/2, (s[:sy] + s[:ey])/2, (s[:sz] + s[:ez])/2

# TODO: Profile this and speed it up: it's probably really slow
"""
    read_sheba(file, origin, stdps=Dict{String,<:Real}) -> s::Vector{SWS}

Read the shear wave splitting measurements `s` from a `file` written by the SHEBA [1]
program.  `origin` is the `Gedoesy.LLA` coordinate of the centre of the UR grid.

To set the depth of stations, supply `stdps` (a `Dict{String,<:Real}`) which
contains the station names as keys and their depths in m as the values.

##### References

1. Wüstefeld, A., Al-Harrasi, O., Verdon, J., Wookey, J. and Kendall, J-M. (2010).
   A strategy for automated analysis of passive microseismic data to image seismic
   anisotropy and fracture characteristics.  Geophysical Prospecting, 58, 755–773.
   doi:10.1111/j.1365-2478.2010.00891.x (available at https://github.com/jwookey/sheba)
"""
function read_sheba(file, origin, stdps=Dict{String,Float64}())
    d = readdlm(file)
    d = d[.!ismatch.(r"^[#%]", string.(d[:,1])), :]
    n = size(d, 1)
    s = Array{SWS}(n)
    for i in 1:n
        elon = d[i,4]
        elat = d[i,3]
        slon = d[i,6]
        slat = d[i,5]
        edep = d[i,7]
        sta = d[i,20]
        phi = d[i,11]
        dphi = d[i,12]
        dt = d[i,13]
        ddt = d[i,14]
        spol = d[i,15]
        dspol = d[i,16]
        stdp = haskey(stdps, sta) ? stdps[sta] : 0.0
        s[i] = SWS(elon, elat, edep, slon, slat, stdp, sta, phi, dphi, dt, ddt,
            spol, dspol, origin)
    end
    s
end

"""
    read_sheba_time(file, args...; comment=r\"^[#%]\") -> st::Vector{SWSTime}

Read results from a SHEBA-formatted file, including the date and time of the
events listed therein, returning the values `st`.

Optionally specify the Regex used to determine which lines are comments with `comment`.
"""
function read_sheba_time(file, args...; comment=r"^[#%]")
    s = read_sheba(file, args...)
    n = length(s)
    st = Vector{SWSTime}(n)
    data = readdlm(file)
    data = data[.!ismatch.(comment, string.(data[:,1])), :]
    size(data, 1) == n || error("Length of splits read ($n) does "*
        "not match the number of non-comment lines ($(size(data, 1))) in file '$file'")
    for i in 1:n
        y = data[i,1]÷1000
        d = data[i,1]%1000
        h = data[i,2]÷100
        m = data[i,2]%100
        date = Date(y) + Day(d)
        st[i] = SWSTime(s[i], DateTime(year(date), month(date), day(date), h, m))
    end
    st
end

"""
    write_sheba(s::Array{<:AbstractSWS}, file)

Write the shear wave splitting measurements `s` to a file `file` in the format used
by the SHEBA [1] program.

##### References

1. Wüstefeld, A., Al-Harrasi, O., Verdon, J., Wookey, J. and Kendall, J-M. (2010).
   A strategy for automated analysis of passive microseismic data to image seismic
   anisotropy and fracture characteristics.  Geophysical Prospecting, 58, 755–773.
   doi:10.1111/j.1365-2478.2010.00891.x (available at https://github.com/jwookey/sheba)
"""
function write_sheba(s::Array{<:AbstractSWS}, file)
    open(file, "w") do f
        println(f, "%  DATE TIME    EVLA    EVLO    STLA    STLO    EVDP    DIST     AZI     BAZ    FAST   DFAST    TLAG   DTLAG    SPOL   DSPOL    WBEG    WEND  STAT %  DATE TIME EIGORIG EIGCORR      Q     SNR   NDF   STAT SNR")
        for ss in s
            datetime = if typeof(ss) <: SWSTime
                t = ss[:time]
                @sprintf("%04d%03d %02d%02d", year(t), _dayofyear(t), hour(t), minute(t))
            else
                "DATE TIME"
            end
            println(f, datetime, " ", ss[:evla], " ", ss[:evlo], " ", ss[:stla], " ",
                ss[:stlo], " ", ss[:evdp], " ",
                SphericalGeom.delta(ss[:evlo], ss[:evla], ss[:stlo], ss[:stla], true), " ",
                azimuth(ss), " ", backazimuth(ss), " ",
                ss[:phi], " ", ss[:dphi], " ", ss[:dt], " ", ss[:ddt], " ", ss[:spol], " ",
                ss[:dspol], " WBEG WEND % ", ss[:stnm], " DATE TIME EIGORIG EIGCORR Q SNR NDF ",
                ss[:stnm], " SNR")
        end
    end
end

"""
    _dayofyear(d::DateTime) -> doy::Int

Return the day of the year of a `DateTime` `d`."""
_dayofyear(d::DateTime) = Dates.value(Date(year(d), month(d), day(d))
                                       - Date(year(d)) + Day(1))

end # module
