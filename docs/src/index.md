# DustExtinction.jl

```@contents
```

## Installation

From the REPL, press `]` to enter Pkg mode
```
(v 1.1) pkg> add DustExtinction
[...]

julia> using DustExtinction
```

## Usage

```jldoctest setup
julia> using DustExtinction

```

### Color laws

```jldoctest setup
julia> ccm89(4000., 3.1)
1.4645557029425842
```

These laws can be applied across higher dimension arrays using the `.` operator

```jldoctest setup
julia> ccm89.([4000., 5000.], 3.1)
2-element Array{Float64,1}:
 1.4645557029425842
 1.122246878899302

```

If you want to apply total extinction $A_V$ it's as simple as multiplcation
```jldoctest setup
julia> a_v=0.3
0.3

julia> a_v * ccm89(4000., 3.1)
0.43936671088277524
```

### Advanced Examples

The color laws also have built-in support for uncertainties using [Measurements.jl](https://github.com/juliaphysics/measurements.jl).

```jldoctest setup
julia> using Measurements

julia> ccm89.([4000. ± 10.5, 5000. ± 10.2], 3.1)
2-element Array{Measurement{Float64},1}:
 1.4646 ± 0.0033
 1.1222 ± 0.003

```

and also support units via [Unitful.jl](https://github.com/painterqubits/unitful.jl) and its subsidiaries. 

```jldoctest setup
julia> using Unitful, UnitfulAstro

julia> mags = ccm89.([4000u"angstrom", 0.5u"μm"], 3.1)
2-element Array{Gain{Unitful.LogInfo{:Magnitude,10,-2.5},:?,Float64},1}:
 1.4645557029425837 mag
  1.122246878899302 mag

julia> fluxes = ones(2)u"erg/s/cm^2"
2-element Array{Quantity{Float64,𝐌*𝐓^-3,Unitful.FreeUnits{(erg, cm^-2, s^-1),𝐌*𝐓^-3,nothing}},1}:
 1.0 erg cm^-2 s^-1
 1.0 erg cm^-2 s^-1

julia> reddened = fluxes .* mags
2-element Array{Quantity{Float64,𝐌*𝐓^-3,Unitful.FreeUnits{(erg, cm^-2, s^-1),𝐌*𝐓^-3,nothing}},1}:
  0.2595241150526157 erg cm^-2 s^-1
 0.35571423768348603 erg cm^-2 s^-1

julia> reddened ≈ @. fluxes * 10 ^ (-0.4 * ustrip(mags))
true

```

You can even combine the two above to get some really nice workflows exploiting all Julia has to offer! This example would shows how you 
could redden some OIR observational data with uncertainties in the flux density.

```jldoctest setup
julia> using Random, Measurements, Unitful, UnitfulAstro

julia> wave = range(0.3, 1.0, length=5)u"μm"
0.3 μm:0.175 μm:1.0 μm

julia> err = randn(MersenneTwister(0), length(wave))u"erg/s/cm^2"
5-element Array{Quantity{Float64,𝐌*𝐓^-3,Unitful.FreeUnits{(erg, cm^-2, s^-1),𝐌*𝐓^-3,nothing}},1}:
   0.6791074260357777 erg cm^-2 s^-1
   0.8284134829000359 erg cm^-2 s^-1
  -0.3530074003005963 erg cm^-2 s^-1
 -0.13485387193052173 erg cm^-2 s^-1
   0.5866170746331097 erg cm^-2 s^-1

julia> flux = @.(ustrip(u"erg/s/cm^2", 300u"erg/s*cm^2" / wave^4) ± ustrip(err))u"erg/s/cm^2"
5-element Array{Quantity{Measurement{Float64},𝐌*𝐓^-3,Unitful.FreeUnits{(erg, cm^-2, s^-1),𝐌*𝐓^-3,nothing}},1}:
  3.7037037037037034e20 ± 0.68 erg cm^-2 s^-1
   5.893140783143162e19 ± 0.83 erg cm^-2 s^-1
 1.6806134238997232e19 ± -0.35 erg cm^-2 s^-1
  6.475979428646598e18 ± -0.13 erg cm^-2 s^-1
                 3.0e18 ± 0.59 erg cm^-2 s^-1

julia> λ_m = @. 0.25 * ccm89(wave, 3.1)
5-element Array{Gain{Unitful.LogInfo{:Magnitude,10,-2.5},:?,Float64},1}:
  0.4545391743613166 mag
  0.3000686820198816 mag
  0.2068677378008488 mag
 0.14069431885378056 mag
 0.10099999999999998 mag

julia> flux .* λ_m
5-element Array{Quantity{Measurement{Float64},𝐌*𝐓^-3,Unitful.FreeUnits{(erg, cm^-2, s^-1),𝐌*𝐓^-3,nothing}},1}:
 2.4368038160281967e20 ± 0.45 erg cm^-2 s^-1
  4.470121662545148e19 ± 0.63 erg cm^-2 s^-1
  1.389059128986988e19 ± 0.29 erg cm^-2 s^-1
  5.688892575540536e18 ± 0.12 erg cm^-2 s^-1
  2.733513699128758e18 ± 0.53 erg cm^-2 s^-1

```

### Dust maps

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


## Reference/API

### Extinction Laws
```@docs
ccm89
od94
cal00
```

### Dust Maps

```@docs
download_sfd98
SFD98Map
ebv_galactic
```
