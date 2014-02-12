# DustExtinction

Models for the interstellar extinction due to dust.

This module exports two functions:

* `ccm89`: Cardelli, Clayton and Mathis (1989) extinction model
* `od94`: O'Donnell (1994) extinction model

## Install

```julia
julia> Pkg.add("git://github.com/kbarbary/DustExtinction.jl.git")
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
* `r_v`: R_V parameter (optional; default is 3.1)

**Returns:**

* `extinction`: Ratio of extinction in magnitudes to A_V
                (extinction at 5494.5 Angstroms)

