# DustExtinction

Models for the interstellar extinction due to dust.

This module exports two functions:

* `ccm89`: [Cardelli, Clayton and Mathis (1989)]
  (http://adsabs.harvard.edu/abs/1989ApJ...345..245C) extinction model.

* `od94`: [O'Donnell (1994)]
  (http://adsabs.harvard.edu/abs/1994ApJ...422..158O)
  extinction model in the optical (between 3030 angstroms and 9090
  angstroms), otherwise identical to Cardelli, Clayton and Mathis (1989).

The implementation of both functions is intended to be as given by the
original papers. An error is raised for values outside the range from
10 to 0.3 inverse microns (1000 to 33333.3 angstroms).

## Install

```julia
julia> Pkg.add("git://github.com/JuliaAstro/DustExtinction.jl.git")
```

## Usage

```julia
julia> using DustExtinction

julia> ccm89(4000., 3.1)
1.4645557029425842

julia> ccm89([4000.; 5000.], 3.1)
2-element Array{Float64,1}:
 1.46456
 1.12225
```

**Parameters:**

* `wave`: wavelength in Angstroms
* `r_v`: R_V parameter

**Returns:**

* `extinction`: Ratio of extinction in magnitudes to A_V
                (extinction in magnitudes at 5494.5 Angstroms)

