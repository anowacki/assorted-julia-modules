__precompile__()

"""
module SphericalGeom provides routines for making calulations on a sphere.

By default, calculations are made in degrees, and all routines have a final
argument `degrees` which may be set to false for radians.  This is advantageous
when making repeated calculations or when using large arrays in terms of speed.
"""
module SphericalGeom

import StaticArrays: SVector

import Base.step

export
    azimuth,
    cart2geog,
    delta,
    geog2cart,
    great_circle,
    great_circle_azimuth,
    nearest_to_gcp,
    project_to_gcp,
    sample,
    step

"""
    cart2geog(x, y, z, degrees::Bool=true) -> lon, lat, r

Compute the longitude, latitude and radius given the cartesian coordinates `x`,
`y` and `z`, where `x` is at (lon,lat) = (0,0), `y` is at (90°,0) and `z` is
through lat = 90°.
"""
function cart2geog(x, y, z, degrees::Bool=true)
    r = sqrt(x^2 + y^2 + z^2)
    r == 0. && return zero(x), zero(x), zero(x)
    lon = atan2(y, x)
    lat = asin(z/r)
    degrees ? (rad2deg(lon), rad2deg(lat), r) : (lon, lat, r)
end

function cart2geog(x::AbstractArray, y::AbstractArray, z::AbstractArray, degrees::Bool=true)
    dims = size(x)
    dims == size(y) == size(z) || throw(ArgumentError("All arrays must have same length"))
    T = promote_type(float.(eltype.((x, y, z)))...)
    lon, lat, r = Array{T}(dims), Array{T}(dims), Array{T}(dims)
    for i in eachindex(lon)
        lon[i], lat[i], r[i] = cart2geog(x[i], y[i], z[i], degrees)
    end
    lon, lat, r
end

"""
    delta(lon1, lat1, lon2, lat2, degrees::Bool=true) -> d

Compute the angular distance `d` on the sphere between two points, (`lon1`, `lat1`)
and (`lon2`, `lat2`).  Points and distance are read and returned in degrees
by default; use `degrees=false` for radians.
"""
function delta(lon1, lat1, lon2, lat2, degrees::Bool=true)
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    d = atan2(sqrt(
               (cos(lat2)*sin(lon2-lon1))^2 + (cos(lat1)*sin(lat2) -
                sin(lat1)*cos(lat2)*cos(lon2-lon1))^2),
               sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1)
              )
    degrees ? rad2deg(d) : d
end

"""
    geog2cart(lon, lat, r, degrees::Bool=true) -> x, y, z
    geog2cart(lon, lat, degrees::Bool=true) -> x, y, z

Return the cartesian coordinates given the geographic longitude, latitude and
radius `lon`, `lat` and `r`.

If `r` is not given, points are returned on the unit sphere.
"""
function geog2cart(lon, lat, r, degrees::Bool=true)
    points_valid(lon, lat, degrees) || error("geog2cart: Points are not on the sphere")
    degrees && begin lon, lat = deg2rad(lon), deg2rad(lat) end
    x = r*cos(lon)*cos(lat)
    y = r*sin(lon)*cos(lat)
    z = r*sin(lat)
    x, y, z
end
function geog2cart(lon::AbstractArray, lat::AbstractArray, r::AbstractArray, degrees::Bool=true)
    size(lon) == size(lat) == size(r) ||
        throw(ArgumentError("Sizes of lon, lat and r must be the same"))
    T = promote_type(float.(eltype.((lon, lat, r)))...)
    x = Array{T}(size(lon))
    y, z = similar(x), similar(x)
    for i in eachindex(lon)
        x[i], y[i], z[i] = geog2cart(lon[i], lat[i], r[i], degrees)
    end
    x, y, z
end
geog2cart(lon::AbstractArray, lat::AbstractArray, degrees::Bool=true) =
    geog2cart(lon, lat, ones(lon), degrees)


"""
    step(lon, lat, az, delta, degrees::Bool=true) -> lon1, lat1

Compute the end point (`lon1`, `lat1`) reached by travelling on the sphere along
azimuth `az` for `delta` angular distance.  Points, angles and distance are read
and returned in degrees by default; use `degrees=false` for radians.
"""
function step(lon, lat, az, delta, degrees::Bool=true)
    if degrees
        delta <= 180 || error("step: delta cannot be more than 180°")
        lon, lat, az, delta = deg2rad(lon), deg2rad(lat), deg2rad(az), deg2rad(delta)
    else
        delta <= pi || error("step: delta cannot be more than π radians")
    end
    lat2 = asin(sin(lat)*cos(delta) + cos(lat)*sin(delta)*cos(az))
    lon2 = lon + atan2(sin(az)*sin(delta)*cos(lat),
                        cos(delta)-sin(lat)*sin(lat2))
    if degrees
        lon2, lat2 = rad2deg(lon2), rad2deg(lat2)
    end
    lon2, lat2
end

"""
    azimuth(lon1, lat1, lon2, lat2, degrees::Bool=true) -> az

Compute the azimuth `az` from point (`lon1`, `lat1`) to (`lon2`, `lat2`) on the sphere.
Points and azimuth are read and returned in degrees by default; use `degrees=false`
for radians.
"""
function azimuth(lon1, lat1, lon2, lat2, degrees::Bool=true)
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    azimuth = atan2(sin(lon2-lon1)*cos(lat2),
                    cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon2-lon1))
    degrees ? rad2deg(azimuth) : azimuth
end

"""
    sample(d=5, degrees::Bool=true) -> lon[], lat[]

Create a set of points which approximately evenly sample the sphere, spaced about `d`
degrees apart (default 5°).  Returns arrays `lon` and `lat` containing the longitude
and latitude (degrees).  Spacing and points are read and returned in degrees by default;
use `false` as last argument to use radians.
"""
function sample(d=5, degrees::Bool=true)
    const nmax = 50000
    lon = Array(Float64, nmax)
    lat = Array(Float64, nmax)
    n = 1
    dlat = d
    dlon = dlat
    lat[n] = 90.0
    lon[n] = 0.0
    lon_i = 0.0
    lat_i = lat[1] - dlat
    while lat_i > -90.0
        dlon_i = dlon/sind(90.0 - lat_i)
        n_i = round(Int, 360.0/dlon_i)
        for i = 1:n_i
            n += 1
            n <= nmax || error("Number of points greater than predetermined limits ($nmax)")
            lat[n] = lat_i
            lon[n] = lon_i
            lon_i = mod(lon_i + dlon_i, 360.0)
        end
        lon_i = mod(lon_i + dlon_i, 360.0)
        lat_i = lat_i - dlat
    end
    n += 1
    lat[n] = -90.0
    lon[n] = 0.0
    lon[1:n] = mod(lon[1:n] + 180.0, 360.0) - 180.0
    degrees ? (lon[1:n], lat[1:n]) : (rad2deg(lon[1:n], rad2deg(lat[1:n])))
end

"""
    project_to_gcp(lon_gc, lat_gc, lon_p, lat_p, degrees::Bool=true) -> lon, lat

Project the point at (`lon_p`, `lat_p`) onto the great circle defined by the
pole to the circle (`lon_gc`, `lat_gc`), returning the point (`lon`, `lat`).
"""
function project_to_gcp(long, latg, lonp, latp, degrees::Bool=true)
    # Convert to vectors
    g = SVector(geog2cart(long, latg, degrees)...)
    p = SVector(geog2cart(lonp, latp, degrees)...)
    # Pole to the gcp containing g and p
    gp = g × p
    # Check for point and pole being the same
    norm(gp) ≈ 0. && error("Point and pole to plane the same")
    gp = normalize(gp)
    pp = gp × g
    acos(p⋅pp) > π/2 && (pp *= -1)
    lon, lat, r = cart2geog(pp[1], pp[2], pp[3], degrees)
    lon, lat
end

"""
    nearest_to_gcp(lon_gc, lat_gc, lons, lats, degrees::Bool=trye) -> i, lon, lat

Find the nearest point to the great circle defined by the pole (`lon_gc`, `lat_gc`)
which occurs in the vectors `lons` and `lats`.  The index `i` and values at that index
`lon` and `lat` are returned.
"""
function nearest_to_gcp(long, latg, lons::AbstractArray, lats::AbstractArray, degrees::Bool=true)
    length(lons) == length(lats) ||
        throw(ArgumentError("lons and lats must have same number of points"))
    mindist = Inf
    imin = 0
    for (i, (lon, lat)) in enumerate(zip(lons, lats))
        lonp, latp = project_to_gcp(long, latg, lon, lat, degrees)
        dist = delta(lon, lat, lonp, latp, degrees)
        if dist < mindist
            mindist = dist
            imin = i
        end
    end
    imin, lons[imin], lats[imin]
end

"""
    great_circle(lon, lat, azimuth, degrees::Bool=true) -> lon_gc, lat_gc
    great_circle(lon1, lat1, lon2, lat2, degrees::Bool=true) -> lon_gc, lat_gc

Determine the great circle, defined by the pole at (`lon_gc`, `lat_gc`), given by
either (the first form), a point and an azimuth, or (the second) two points on the
sphere.
"""
function great_circle(lon, lat, azimuth, degrees::Bool=true)
    dist = degrees ? 45 : π/4
    lonp, latp = step(lon, lat, azimuth, dist, degrees)
    great_circle(lon, lat, lonp, latp, degrees)
end
function great_circle(lon1, lat1, lon2, lat2, degrees::Bool=true)
    p1 = SVector(geog2cart(lon1, lat1, degrees)...)
    p2 = SVector(geog2cart(lon2, lat2, degrees)...)
    p1 ≈ p2 && error("Two points overlap (are: $p1 and $p2)")
    delta(lon1, lat1, lon2, lat2, degrees) ≈ (degrees ? 180 : π) &&
        error("Points are antipodal (are: $p1 and $p2)")
    gp = p1 × p2
    lon, lat, r = cart2geog(gp[1], gp[2], gp[3], degrees)
    lon, lat
end

"""
    great_circle_azimuth(long, latg, lonp, latp, degrees::Bool=true) -> azimuth

Return the local `azimuth` of the great circle with pole (`long`, `latg`) where
the point (`lonp`, `latp`) projects onto the great circle.
"""
function great_circle_azimuth(long, latg, lonp, latp, degrees::Bool=true)
    lon, lat = project_to_gcp(long, latg, lonp, latp, degrees)
    mod(azimuth(lon, lat, long, latg, degrees) + (degrees ? 90 : π/2), (degrees ? 360 : 2π))
end

"""
    points_valid(lon, lat, degrees::Bool=true) -> ::Bool

Return `true` if all points in arrays `lon` and `lat` are on the sphere.
`lon` is ignored, but `lat` is checked to see if points are in the range
-90°–90° (-π–π).  Points are read and returned in degrees by default; use
`degrees=false` for radians.
"""
points_valid(lon, lat, degrees::Bool=true) =
    degrees ? !any(abs.(lat) .> 90.) : !any(abs.(lat) .> pi/2.)

end # module
