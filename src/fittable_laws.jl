

"""
    FM90(;c1=0.10, c2=0.70, c3=3.23, c4=0.41, x0=4.60, gamma=0.9)(λ::Union{Real, Quantity})
    FM90(coeffs, x0=4.60, gamma=0.9)(λ::Union{Real, Quantity})

Fitzpatrick & Massa (1990) 6 parameter ultraviolet shape model, this model is only applicable at UV wavelengths.

### Parameters
* `c1` - y-intercept of linear term
* `c2` - slope of liner term
* `c3` - amplitude of 2175 Å bump
* `c4` - amplitude of FUV rise
* `x0` - centroid of 2175 Å bump
* `gamma` - width of 2175 Å bump

The model has form ``c_1 + c_2x + c_3D(x; \\gamma, x_0) + c_4 F(x)`` where ``x`` is the wavenumber in inverse microns, ``D(x)``
is a Drude profile (modified Lorentzian) used to model the 2175 Å bump with the scale-free parameters ``x_0`` (central wavenumber)
and ``\\gamma`` (damping coefficient), and ``F(x)``, a piecewise function for the far-UV. Note that the coefficients will change
the overall normalization, possibly changing the expected behavior of reddening via the parameter ``A_V``.

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the
returned value will be `UnitfulAstro.mag`.

# Examples
```jldoctest
julia> model = FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99);

julia> model(1500)
5.2521258452800135

julia> FM90()(1500)
5.152125845280013

julia> FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99).([1000, 1200, 1800])
3-element Array{Float64,1}:
 12.562237969522851
  7.769215017329513
  4.890128210972148

```
"""
@with_kw struct FM90{T<:Number} @deftype T
    c1 = 0.10
    c2 = 0.70
    c3 = 3.23
    c4 = 0.41
    x0 = 4.60
    gamma = 0.99
    @assert x0 ≥ 0 "`x0` must be ≥ 0, got $x0"
    @assert gamma ≥ 0 "`gamma` must be ≥ 0, got $gamma"
end

FM90(c1, c2, c3, c4, x0, gamma) = FM90(promote(c1, c2, c3, c4, x0, gamma)...)
FM90(coeffs, x0=4.60, gamma=0.99) = FM90(coeffs..., x0, gamma)

function (z::FM90)(λ::Real)
    x = aa_to_invum(λ)

    if !(1.0 / 0.32 <= x <= 1.0 / 0.0912)
        return zero(float(eltype(x)))
    end

    exvebv = z.c1 + z.c2 * x

    x2 = x^2
    exvebv += z.c3 * (x2 / ((x2 - z.x0^2)^2 + x2 * (z.gamma^2)))

    if x >= 5.9
        y = x - 5.9
        exvebv += z.c4 * @evalpoly y 0.0 0.0 0.5392 0.05644
    end

    return exvebv
end

(z::FM90)(λ::Quantity) = z(ustrip(u"angstrom", λ)) * u"mag"
