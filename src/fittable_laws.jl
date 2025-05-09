
"""
    FM90(;c1=0.10, c2=0.70, c3=3.23, c4=0.41, x0=4.60, gamma=0.9)
    FM90(coeffs, x0=4.60, gamma=0.9)

Fitzpatrick & Massa (1990) generative model for ultraviolet dust extinction. The
default values are the published values for the Milky Way average.

## Parameters
* `c1` - y-intercept of linear term
* `c2` - slope of liner term
* `c3` - amplitude of 2175 Å bump
* `c4` - amplitude of FUV rise
* `x0` - centroid of 2175 Å bump
* `gamma` - width of 2175 Å bump

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the
returned value will be `UnitfulAstro.mag`.

## Examples
```jldoctest
julia> model = FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99);

julia> model(1500)
5.2521258452800135

julia> FM90()(1500)
5.152125845280013

julia> FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99).([1000, 1200, 1800])
3-element Vector{Float64}:
 12.562237969522851
  7.769215017329513
  4.890128210972148

```
# Extended Help

The model has form ``c_1 + c_2x + c_3D(x; γ, x_0) + c_4 F(x)`` where ``x`` is
the wavenumber in inverse microns, ``D(x)`` is a Drude profile (modified Lorentzian)
used to model the 2175 Å bump with the scale-free parameters ``x_0`` (central wavenumber)
and ``γ`` (damping coefficient), and ``F(x)``, a piecewise function for the far-UV.
Note that the coefficients will change the overall normalization, possibly
changing the expected behavior of reddening via the parameter ``A_V``.

## References
[Fitzpatrick & Massa (1990)](https://ui.adsabs.harvard.edu/abs/1990ApJS...72..163F)
"""
Base.@kwdef struct FM90{T<:Number} <: ExtinctionLaw
    c1::T = 0.10
    c2::T = 0.70
    c3::T = 3.23
    c4::T = 0.41
    x0::T = 4.60
    gamma::T = 0.99
    function FM90(c1, c2, c3, c4, x0, gamma)
        x0 < 0 && error("`x0` must be ≥ 0, got ", x0)
        gamma < 0 && error("`gamma` must be ≥ 0, got ", gamma)
        params = promote(c1, c2, c3, c4, x0, gamma)
        return new{eltype(params)}(params...)
    end
end
FM90(coeffs, x0=4.60, gamma=0.99) = FM90(coeffs..., x0, gamma)

bounds(::Type{<:FM90}) = (912, 3200)

function (law::FM90)(wave::T) where T <: Real
    checkbounds(law, wave) || return zero(float(T))

    x = aa_to_invum(wave)
    exvebv = law.c1 + law.c2 * x

    x2 = x^2
    exvebv += law.c3 * (x2 / ((x2 - law.x0^2)^2 + x2 * (law.gamma^2)))

    if x >= 5.9
        y = x - 5.9
        exvebv += law.c4 * @evalpoly y 0 0 0.5392 0.05644
    end

    return exvebv
end
