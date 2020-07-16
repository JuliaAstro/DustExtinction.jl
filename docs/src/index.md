# DustExtinction.jl

This package implements different empirical dust measurements for use in astronomy. This package is written in pure Julia and is built with first-class support with [Unitful.jl](https://github.com/painterqubits/unitful.jl) and [Measurements.jl](https://github.com/juliaphysics/measurements.jl).

## About
Extinction describes the effect of dust grains on observations of stars in space. Light that travels through dust is absorbed and scatterred as natural processes of light's interactions with materials. This obfuscation can be modeled and removed from our data in order to more properly retrieve the star's flux. When dealing with multiple stars, or large clusters or galaxies, this process is considered dust attenuation and is not provided for in this package.

## Installation

From the REPL, press `]` to enter Pkg mode

```julia
(v 1.4) pkg> add DustExtinction

julia> using DustExtinction
```

## Usage

```jldoctest
julia> using DustExtinction

julia> CCM89(Rv=3.1)(4000)
1.464555702942584
```

For more examples, view the [Color Laws](@ref laws) and [Dust Maps](@ref maps) sections.

## Citations

There are various citations relevant to this work. Please be considerate when using this work or any derivate of it by adding the appropriate citations.

|        Law         | Reference                                                                                      | BibTeX                       |
| :----------------: | :--------------------------------------------------------------------------------------------- | :--------------------------- |
|  [`CCM89`](@ref)   | [Clayton, Cardelli and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C)   | [download](assets/ccm89.bib) |
|   [`OD94`](@ref)   | [O'Donnell (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O)                      | [download](assets/od94.bib)  |
|  [`CAL00`](@ref)   | [Calzetti et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C)                | [download](assets/cal00.bib) |
|  [`VCG04`](@ref)   | [Valencic, Clayton, & Gordon (2004)](https://ui.adsabs.harvard.edu/abs/2004ApJ...616..912V)    | [download](assets/vcg04.bib) |
|  [`GCC09`](@ref)   | [Gordon, Cartledge, & Clayton (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...705.1320G)   | [download](assets/gcc09.bib) |
|  [`FM90`](@ref)   | [Fitzpatrick & Massa (1990)](https://ui.adsabs.harvard.edu/abs/1990ApJS...72..163F)   | [download](assets/fm90.bib) |
| [`SFD98Map`](@ref) | [Schlegel, Finkbeiner and Davis (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S) | [download](assets/sfd98.bib) |
| [`F99`](@ref) | [Fitzpatrick (1999)](https://ui.adsabs.harvard.edu/abs/1999PASP..111...63F) | [download](assets/f99.bib) |
| [`M14`](@ref) | [Maiz Apellaniz et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014A%26A...564A..63M) | [download](assets/m14.bib) |

## Index

```@index
Modules = [DustExtinction]
Order = [:function, :type]
```

## Contributing

If you are interested in contributing, feel free to make a pull request or open an issue for discussion.
