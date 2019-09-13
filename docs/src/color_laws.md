
# [Color laws](@id laws)

The following empirical laws allow us to model the reddening of light as it travels to us. If you aren't quite sure which law to use, [`ccm89`](@ref) is a decent default.

```@meta
DocTestSetup = quote
    using DustExtinction, Random
    Random.seed!(1)
end
```

```@raw html
<div style="">
<img src="assets/laws.png" style="max-width:80%; display: block; margin-left: auto; margin-right: auto">
</div>
```
Here is a comparison of these three laws in the 0.3-10 inverse-micron range. Notice the slight difference in `od94` and `ccm89` around the optical (2.5-3 inverse-micron) while they are equal everywhere else.

## Usage

```jldoctest
julia> ccm89(4000., 3.1)
1.4645557029425842

```

These laws can be applied across higher dimension arrays using the `.` operator

```jldoctest
julia> ccm89.([4000., 5000.], 3.1)
2-element Array{Float64,1}:
 1.4645557029425842
 1.122246878899302

```

these laws return magnitudes, which we can apply directly to flux by mulitplication with a base-2.5 logarithmic system (because astronomers are fun):

```julia
f = f * 10 ^ (-0.4 * Av * mag)
```

To make this easier, we provide a convenient [`extinct`](@ref) function for applying these color laws to flux measurements.

```jldoctest
julia> wave = range(4000, 5000, length=4)
4000.0:333.3333333333333:5000.0

julia> flux = 1e-8 .* wave .+ 1e-2
0.01004:3.3333333333333333e-6:0.01005

julia> extinct.(flux, wave, 0.3)
4-element Array{Float64,1}:
 0.006698646015454752 
 0.006918253926353551 
 0.007154659823737299 
 0.0073704912727315395

```

!!! note
    If you are concerned with memory and would like an in-place broadcasted version of [`extinct`](@ref) use it with [`Base.map!`](https://docs.julialang.org/en/v1/base/collections/#Base.map!)
    ```julia
    julia> map!((f,w)->extinct(f, w, 0.3), flux, flux, wave)
    ```

## Advanced Usage

The color laws also have built-in support for uncertainties using [Measurements.jl](https://github.com/juliaphysics/measurements.jl).

```jldoctest setup
julia> using Measurements

julia> ccm89.([4000. ± 10.5, 5000. ± 10.2], 3.1)
2-element Array{Measurement{Float64},1}:
 1.4646 ± 0.0033
 1.1222 ± 0.003

```

and also support units via [Unitful.jl](https://github.com/painterqubits/unitful.jl) and its subsidiaries. Notice how the output type is now a `Unitful.Gain`.

```jldoctest
julia> using Unitful, UnitfulAstro

julia> mags = ccm89.([4000u"angstrom", 0.5u"μm"], 3.1)
2-element Array{Gain{Unitful.LogInfo{:Magnitude,10,-2.5},:?,Float64},1}:
 1.4645557029425837 mag
  1.122246878899302 mag

```

You can even combine the two above to get some really nice workflows exploiting all Julia has to offer! This example shows how you 
could redden some OIR observational data with uncertainties in the flux density.

```jldoctest setup
julia> using Measurements, Unitful, UnitfulAstro

julia> wave = range(0.3, 1.0, length=5)u"μm"
(0.3:0.175:1.0) μm

julia> err = randn(length(wave))
5-element Array{Float64,1}:
  0.2972879845354616
  0.3823959677906078
 -0.5976344767282311
 -0.01044524463737564
 -0.839026854388764

julia> flux = @.(300 / ustrip(wave)^4 ± err)*u"Jy"
5-element Array{Quantity{Measurement{Float64},𝐌*𝐓^-2,Unitful.FreeUnits{(Jy,),𝐌*𝐓^-2,nothing}},1}:
  37037.04 ± 0.3 Jy
  5893.14 ± 0.38 Jy
  1680.61 ± -0.6 Jy
 647.598 ± -0.01 Jy
   300.0 ± -0.84 Jy

julia> extinct.(flux, wave, 0.3)
5-element Array{Quantity{Measurement{Float64},𝐌*𝐓^-2,Unitful.FreeUnits{(Jy,),𝐌*𝐓^-2,nothing}},1}:
    22410.8 ± 0.18 Jy
    4229.74 ± 0.27 Jy
    1337.12 ± 0.48 Jy
 554.3349 ± 0.0089 Jy
     268.31 ± 0.75 Jy

```

## API/Reference

```@docs
extinct
ccm89
od94
cal00
```
