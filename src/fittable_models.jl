# Convenience function for wavelength conversion
aa_to_invum(wave::Real) = 10000 / wave
aa_to_invum(wave::Quantity) = aa_to_invum(ustrip(u"angstrom", wave))

"""
    FM90(C1 = 0.10, C2 = 0.70, C3 = 3.23, C4 = 0.41, xo = 4.60, γ = 0.9)(λ::Real)
    FM90(C1 = 0.10, C2 = 0.70, C3 = 3.23, C4 = 0.41, xo = 4.60, γ = 0.9)(λ::Quantity)

### Parameters
* `C1` - y-intercept of linear term
* `C2` - slope of liner term
* `C3` - amplitude of 2175 Å bump
* `C4` - amplitude of FUV rise
* `xo` - centroid of 2175 Å bump
*  `γ` - width of 2175 Å bump

Fitzpatrick & Massa (1990) 6 parameter ultraviolet shape model, this model is only applicable at UV wavelengths.

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the
returned value will be `UnitfulAstro.mag`.
"""
struct FM90{T <: Number}
    C1::T
    C2::T
    C3::T
    C4::T
    xo::T
    γ::T

    function FM90(C1::T, C2::T, C3::T, C4::T, xo::T, γ::T) where T <: Number
        xo < 0 && error("Invalid axis xo=$xo. xo must be greater than or equal to 0")
        γ < 0 && error("Invalid axis γ=$γ. γ must be greater than or equal to 0")
        new{T}(C1, C2, C3, C4, xo, γ)
    end
end

FM90(C1, C2, C3, C4, xo, γ) = FM90(promote(C1, C2, C3, C4, xo, γ)...)
FM90() = FM90(0.10, 0.70, 3.23, 0.41, 4.60, 0.99)

function (z::FM90)(λ::Real)
    x = aa_to_invum(λ)

    if !(1.0 / 0.32 <= x <= 1.0 / 0.0912)
        return zero(float(eltype(x)))
    end

    exvebv = z.C1 + z.C2 * x

    x2 = x^2
    exvebv += z.C3 * (x2 / ((x2 - z.xo^2)^2 + x2 * (z.γ^2)))

    if x >= 5.9
        y = x - 5.9
        exvebv += z.C4 * @evalpoly y 0.0 0.0 0.5392 0.05644
    end

    return exvebv
end

(z::FM90)(λ::Quantity) = z(ustrip(u"angstrom", λ)) * u"mag"
