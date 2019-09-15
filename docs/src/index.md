# DustExtinction.jl

This package implements different empirical dust measurements for use in astronomy. This package is written in pure Julia and is built with first-class support with [`Unitful.jl`](https://github.com/painterqubits/unitful.jl) and [`Measurements.jl`](https://github.com/juliaphysics/measurements.jl).

## About
Extinction describes the effect of dust grains on observations of stars in space. Light that travels through dust is absorbed and scatterred as natural processes of light's interactions with materials. This obfuscation can be modeled and removed from our data in order to more properly retrieve the star's flux. When dealing with multiple stars, or large clusters or galaxies, this process is considered dust attenuation and is not provided for in this package.

## Installation

From the REPL, press `]` to enter Pkg mode

```julia
(v 1.1) pkg> add DustExtinction
[...]

julia> using DustExtinction
```

## Usage

```jldoctest
julia> using DustExtinction

julia> ccm89(4000., 3.1)
ERROR: MethodError: no method matching ccm89_invum(::Float64, ::Float64, ::Array{Float64,1}, ::Array{Float64,1})
Closest candidates are:
  ccm89_invum(::Real, ::Real, !Matched::Polynomials.Poly{#s17} where #s17<:Real, !Matched::Polynomials.Poly{#s16} where #s16<:Real) at /Users/miles/dev/julia/DustExtinction.jl/src/ccm89.jl:10
Stacktrace:
 [1] ccm89(::Float64, ::Float64) at /Users/miles/dev/julia/DustExtinction.jl/src/ccm89.jl:59
 [2] top-level scope at none:4

```

For more examples, view the [Color Laws](@ref laws) and [Dust Maps](@ref maps) sections.

## Citations

There are various citations relevant to this work. Please be considerate when using this work or any derivate of it by adding the appropriate citations.

[`ccm89`](@ref) - 
[Cardelli, Clayton and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract) - [BibTex](assets/ccm89.bib)

[`od94`](@ref) - [O'Donnell (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O/abstract) - [BibTex](assets/od94.bib)

[`cal00`](@ref) - [Calzetti et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract) - [BibTex](assets/cal00.bib)

[`SFD98Map`](@ref) - [Schlegel, Finkbeiner and Davis (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract) - [BibTex](assets/sfd98.bib)

## Index

```@index
Modules = [DustExtinction]
Order = [:function, :type]
```

## Contributing

If you are interested in contributing, feel free to make a pull request or open an issue for discussion.
