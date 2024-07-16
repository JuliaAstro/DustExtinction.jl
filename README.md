# DustExtinction.jl


[![Build Status](https://github.com/juliaastro/DustExtinction.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/DustExtinction.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/D/DustExtinction.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/DustExtinction.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/DustExtinction.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/DustExtinction.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/DustExtinction.jl/dev)



Tools for interstellar dust extinction in astronomy

Documentation: [DustExtinction](https://JuliaAstro.github.io/DustExtinction.jl/stable)

## Installation

From package manager (`]` key)

```julia-repl
pkg> add DustExtinction
```

## Usage

Color laws:

```julia-repl
julia> using DustExtinction

julia> CCM89(Rv=3.1)(4000)
1.4645557029425842

julia> CCM89(Rv=3.1).([4000, 5000])
2-element Vector{Float64}:
 1.46456
 1.12225
```

Dust maps:

```julia-repl
julia> dustmap = SFD98Map()
SFD98Map("[...]")

julia> dustmap(0.1, 0.1)
0.793093095733043

julia> dustmap.([0.1, 0.2], [0.1, 0.2])
2-element Vector{Float64}:
 0.793093
 0.539507
```

Reddening:

```julia-repl
julia> wave = [4000., 5000.]
2-element Vector{Float64}:
 4000.0
 5000.0

julia> flux = [1.0, 1.5]
2-element Vector{Float64}:
 1.0
 1.5

julia> red = redden.(CCM89, wave, flux; Av=0.3, Rv=3.1)
2-element Vector{Float64}:
 0.6671958182723856
 1.1000733242882896

julia> deredden.(CCM89(Rv=3.1), wave, red; Av=0.3)
2-element Vector{Float64}:
 1.0
 1.5
```

We provide first-class support for `Unitful.jl` and `Measurements.jl` packages, too! Check out the documentation for more examples.

## Contributing

Feel free to open an issue or a pull-request for any discussion, suggestions, new features, or patches!
