# DustExtinction.jl

[![Build Status](https://img.shields.io/travis/JuliaAstro/DustExtinction.jl.svg)](https://travis-ci.org/JuliaAstro/DustExtinction.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/DustExtinction.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/DustExtinction.jl/dev)
[![codecov](https://codecov.io/gh/JuliaAstro/DustExtinction.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/DustExtinction.jl)



Tools for interstellar dust extinction in astronomy

Documentation: [DustExtinction](https://JuliaAstro.github.io/DustExtinction.jl/stable)

## Installation

From package manager (``` ]``` key)

``` add DustExtinction```

## Usage

Color laws:

```julia
julia> using DustExtinction

julia> ccm89(4000., 3.1)
1.4645557029425842

julia> ccm89.([4000., 5000.], 3.1)
2-element Array{Float64,1}:
 1.46456
 1.12225
```

Dust maps:

```julia
julia> dustmap = SFD98Map()
SFD98Map("[...]")

julia> dustmap(0.1, 0.1)
0.793093095733043

julia> dustmap.([0.1, 0.2], [0.1, 0.2])
2-element Array{Float64,1}:
 0.793093
 0.539507
```

Extinction:

```julia
julia> wave = [4000., 5000.]
2-element Array{Float64,1}:
 4000.0
 5000.0

julia> flux = [1.0, 1.5]
2-element Array{Float64,1}:
 1.0
 1.5

julia> extinct.(flux, wave, 0.3)
2-element Array{Float64,1}:
 0.6671958182723856
 1.1000733242882896

```

We provide first-class support for `Unitful.jl` and `Measurements.jl` packages, too! Check out the documentation for more examples.

## Contributing

Feel free to open an issue or a pull-request for any discussion, suggestions, new features, or patches!
