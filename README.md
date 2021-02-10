# `assorted-julia-modules`

An old repo containing a few useful Julia modules.

## Abandonment warning
**Most of the modules in this repo are abandoned or deprecated.
Some have been replaced by packages available at
https://github.com/anowacki.  Some still work.
None of them come with any
guarantee of correctness or usability, nor any promise to fix
bugs.**

## Usable modules
These modules have been updated for Julia v1 and so should still work.

- `CircPlot`: Plot circular data, e.g. polar histograms.

  [**Note:**
  This is the only package I would be prepared to help with, say
  if someone wants to add these plots to something like
  [StatsPlots](https://github.com/JuliaPlots/StatsPlots.jl).
  ]

- `CorrelationDimension`: Calculate the correlation or fractal dimension
  of a set of points
- `FortranReader`: Read and write 'non-portable' Fortran unformatted files
  produced by at least `gfortran`.
- `SphericalGeom`: Routines for points on the sphere.

## Superceded modules

The following modules have been directly replaced by properly
tested packages:
- `AK135`, `PREM`: Use [SeisModels](https://github.com/anowacki/SeisModels.jl).

## Installing

Clone this repository somewhere, then add the directory to your
`LOAD_PATH` by doing something like
```julia
shell> git clone https://github.com/anowacki/assorted-julia-modules

julia> push!(LOAD_PATH, joinpath(pwd(), "assorted-julia-modules"))
```
Note that this is **not** the recommended way to add
Julia packages to your environment.

You will then need to manually add the dependencies for each package:

- `CircPlot`: `pkg> add StatsBase Plots`
- `CorrelationDimension`: `pkg> add Plots Distributions StatsBase`
- `SphericalGeom`: `pkg> add StaticArrays`

## License
See [the licence file](LICENCE.md).
