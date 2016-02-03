"""
module Iso3D extends PyPlot to easily create isosurface plots from regularly-
spaced 3d arrays of data.
It requires matplotlib, numpy and scikit-image be importable into the default Python.
"""
module Iso3D

using PyCall, PyPlot

export isosurface

@pyimport mpl_toolkits.mplot3d.art3d as art3d
@pyimport skimage.measure as measure

"""
`isosurface(a::Array{T,3}, v)`

Draw a 3D isosurface at value `v` drawn from the 3D array `a`.
"""
function isosurface{T}(a::Array{T,3}, v)
	minimum(a) < v && maximum(a) > v || error("Isosurface value is not contained in data")
	verts, faces = measure.marching_cubes(a, v)
	nf = size(faces, 1)
	triangles = Array{eltype(verts)}(nf, 3, 3)
	x = verts[:,1]'
	y = verts[:,2]'
	z = verts[:,3]'
	@inbounds for i = 1:nf, j = 1:3, k = 1:3
		triangles[i,j,k] = verts[faces[i,j]+1,k]
	end
	collection = art3d.Poly3DCollection(triangles)
	fig = figure(figsize=(10, 10))
	ax = subplot(111, projection="3d")
	ax[:add_collection](collection)
	xlim(minimum(x), maximum(x))
	ylim(minimum(y), maximum(y))
	zlim(minimum(z), maximum(z))
	xlabel("x")
	ylabel("y")
	zlabel("z")
end

end
