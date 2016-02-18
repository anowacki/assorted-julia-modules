"""
module SphericalGeom provides routines for making calulations on a sphere.

By default, calculations are made in degrees, and all routines have a final
argument `degrees` which may be set to false for radians.  This is advantageous
when making repeated calculations or when using large arrays in terms of speed.
"""
module SphericalGeom

export
    azimuth,
    delta,
    step

"""
`delta(lon1, lat1, lon2, lat2, degrees::Bool=true) -> d`

Compute the angular distance `d` on the sphere between two points, (`lon1`, `lat1`)
and (`lon2`, `lat2`).  Points and distance are read and returned in degrees
by default; use `degrees=false` for radians.
"""
function delta(lon1, lat1, lon2, lat2, degrees::Bool=true)
    points_valid([lon1, lon2], [lat1, lat2], degrees) ||
        error("delta: Points are not on the sphere")
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
`step(lon, lat, az, delta, degrees::Bool=true) -> lon1, lat1`

Compute the end point (`lon1`, `lat1`) reached by travelling on the sphere along
azimuth `az` for `delta` angular distance.  Points, angles and distance are read
and returned in degrees by default; use `degrees=false` for radians.
"""
function step(lon, lat, az, delta, degrees::Bool=true)
    points_valid(lon, lat, degrees) ||
        error("step: Points are not on the sphere")
    if degrees
        all(delta .<= 180.) || error("step: delta cannot be more than 180°")
        lon, lat, az, delta = deg2rad(lon), deg2rad(lat), deg2rad(az), deg2rad(delta)
    else
        all(delta .<= pi) || error("step: delta cannot be more than π radians")
    end
    lat2 = asin(sin(lat).*cos(delta) + cos(lat).*sin(delta).*cos(az))
    lon2 = lon + atan2(sin(az).*sin(delta).*cos(lat),
                        cos(delta)-sin(lat).*sin(lat2))
    if degrees
        lon2, lat2 = rad2deg(lon2), rad2deg(lat2)
    end
    lon2, lat2
end

"""
`azimuth(lon1, lat1, lon2, lat2, degrees::Bool=true) -> az`

Compute the azimuth `az` from point (`lon1`, `lat1`) to (`lon2`, `lat2`) on the sphere.
Points and azimuth are read and returned in degrees by default; use `degrees=false`
for radians.
"""
function azimuth(lon1, lat1, lon2, lat2, degrees::Bool=true)
    points_valid([lon1, lon2], [lat1, lat2], degrees) ||
        error("azimuth: Points are not on the sphere")
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    azimuth = atan2(sin(lon2-lon1).*cos(lat2),
                    cos(lat1).*sin(lat2) - sin(lat1).*cos(lat2).*cos(lon2-lon1))
    degrees ? rad2deg(azimuth) : azimuth
end

"""
`sample(d=5, degrees::Bool=true) -> lon[], lat[]`

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
    lat[n] = 90.
    lon[n] = 0.
    lon_i = 0.
    lat_i = lat[1] - dlat
    while lat_i > -90.
        dlon_i = dlon/sind(90. - lat_i)
        n_i = round(Int, 360./dlon_i)
        for i = 1:n_i
            n += 1
            n <= nmax || error("Number of points greater than predetermined limits ($nmax)")
            lat[n] = lat_i
            lon[n] = lon_i
            lon_i = mod(lon_i + dlon_i, 360.)
        end
        lon_i = mod(lon_i + dlon_i, 360.)
        lat_i = lat_i - dlat
    end
    n += 1
    lat[n] = -90.
    lon[n] = 0.
    lon[1:n] = mod(lon[1:n] + 180., 360.) - 180.
    degrees ? (lon[1:n], lat[1:n]) : (rad2deg(lon[1:n], rad2deg(lat[1:n])))
end


"""
`points_valid(lon, lat, degrees::Bool=true) -> ::Bool`

Return `true` if all points in arrays `lon` and `lat` are on the sphere.
`lon` is ignored, but `lat` is checked to see if points are in the range
-90°–90° (-π–π).  Points are read and returned in degrees by default; use
`degrees=false` for radians.
"""
points_valid(lon, lat, degrees::Bool=true) =
    degrees ? all(abs(lat) .<= 90.) : all(abs(lat) .<= pi)

end # module

