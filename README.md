# DustExtinction.jl

[![Build Status](https://img.shields.io/travis/JuliaAstro/DustExtinction.jl.svg?style=flat-square)](https://travis-ci.org/JuliaAstro/DustExtinction.jl)


Tools for interstellar dust extinction in astronomy

## Install

From package manager (``` ]``` key)

``` add DustExtinction```

Usage
-----

Color laws:

```julia
julia> using DustExtinction

julia> ccm89(4000., 3.1)
1.4645557029425842

julia> ccm89([4000., 5000.], 3.1)
2-element Array{Float64,1}:
 1.46456
 1.12225
```

Dust maps:

```julia
julia> ENV["SFD98_DIR"] = "/home/user/data/dust"

# download maps (once)
julia> download_sfd98()

julia> dustmap = SFD98Map()
SFD98Map("/home/user/data/dust")

julia> ebv_galactic(dustmap, 0.1, 0.1)
0.793093095733043

julia> ebv_galactic(dustmap, [0.1, 0.2], [0.1, 0.2])
2-element Array{Float64,1}:
 0.793093
 0.539507
```


Reference/API
-------------

* `ccm89(wave, r_v)`

  [Cardelli, Clayton and Mathis (1989)]
  (http://adsabs.harvard.edu/abs/1989ApJ...345..245C)
  extinction function. Given wavelength(s) in Angstroms, returns
  A(wave)/A_V: extinction in magnitudes relative to value at 5494.5
  Angstroms.

  The implementation of this function and `od94` is intended to be as
  given by the original papers. An error is raised for values outside
  the range from 10 to 0.3 inverse microns (1000 to 33333.3
  angstroms).

* `od94(wave, r_v)`

  [O'Donnell (1994)](http://adsabs.harvard.edu/abs/1994ApJ...422..158O)
  extinction model in the optical (between 3030 angstroms and 9090
  angstroms), otherwise identical to Cardelli, Clayton and Mathis (1989).

* `download_sfd98([destdir])`

  Download the Schlegel, Finkbeiner and Davis (1998) galactic dust map to
  the local directory `destdir`, if specified. The default destination
  directory is the `$SFD98_DIR` environment variable.

* `SFD98Map([mapdir])`

  Type to represent the SFD (1998) map. If `mapdir` not given, default is
  the `$SFD98_DIR` variable.

* `ebv_galactic(dustmap, l, b)`

  Retrieve *E(B-V)* value from an SFD98Map at given galactic longitude and
  latitude. Return value is linearly interpolated between pixel positions in
  the map.
