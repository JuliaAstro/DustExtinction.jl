# Convenience function for wavelength conversion
aa_to_invum(wave::Real) = 10000 / wave
aa_to_invum(wave::Quantity) = aa_to_invum(ustrip(u"angstrom", wave))

"""
    FM90(C1 = 0.10, C2 = 0.70, C3 = 3.23, C4 = 0.41, xo = 4.60, γ = 0.9)(λ::Real)

The parameters C1, C2, C3, C4, xo and γ are y-intercept, slope, bump amplitude, FUV rise amplitude,
bump centroid and bump width respectively.

Fitzpatrick & Massa (1990) 6 parameter ultraviolet shape model, this model is only applicable at UV wavelengths.

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the
returned value will be `UnitfulAstro.mag`.
"""
struct FM90
    C1::Real
    C2::Real
    C3::Real
    C4::Real
    xo::Real
    γ::Real
end

FM90() = FM90(0.10, 0.70, 3.23, 0.41, 4.60, 0.9)

function (z::FM90)(λ::Real)
    x = aa_to_invum(λ)

    if !(1.0 / 0.32 <= x <= 1.0 / 0.0912)
        return zero(eltype(x))
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
