# Convenience function for wavelength conversion
aa_to_invum(wave::Real) = 10000 / wave
aa_to_invum(wave::Quantity) = aa_to_invum(ustrip(u"angstrom", wave))

"""
    FM90(λ::Real, C1 = 0.10, C2 = 0.70, C3 = 3.23, C4 = 0.41, xo = 4.60, γ = 0.9)
    FM90(λ::Quantity, C1 = 0.10, C2 = 0.70, C3 = 3.23, C4 = 0.41, xo = 4.60, γ = 0.9)

The parameters C1, C2, C3, C4, xo and gamma are y-intercept, slope, bump amplitude, FUV rise amplitude,
bump centroid and bump width respectively.

Fitzpatrick & Massa (1990) 6 parameter ultraviolet shape model, this model is only applicable at UV wavelengths.

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the
returned value will be `UnitfulAstro.mag`.
"""
function FM90(λ::Real, C1 = 0.10, C2 = 0.70, C3 = 3.23, C4 = 0.41, xo = 4.60, γ = 0.9)
    x = aa_to_invum(λ)

    if !(1.0 / 0.32 <= x <= 1.0 / 0.0912)
        return zero(eltype(x))
    end

    exvebv = C1 + C2 * x

    x2 = x^2
    exvebv += C3 * (x2 / ((x2 - xo^2)^2 + x2 * (γ^2)))

    if x >= 5.9
        y = x - 5.9
        exvebv += C4 * @evalpoly y 0.0 0.0 0.5392 0.05644
    end

    return exvebv
end

FM90(λ::Quantity, C1::Real = 0.10,
     C2::Real = 0.70, C3::Real = 3.23,
     C4::Real = 0.41, xo::Real = 4.60,
     γ::Real = 0.9) = FM90(ustrip(u"angstrom", λ), C1, C2, C3, C4, xo, γ) * u"mag"
