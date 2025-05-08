# Optical
const m14_x1 = (1.0)
const m14_xi1 = first(m14_x1)
const m14_x2 = (1.15, 1.81984, 2.1, 2.27015, 2.7)
const m14_x3 = (3.5, 3.9, 4.0, 4.1, 4.2)
const m14_xi3 = last(m14_x3)
const m14_a1v = 0.574 * m14_x1^1.61
const m14_b1v = -0.527 * m14_x1^1.61
const m14_x2_diff = m14_x2 .- 1.82
const m14_a2v = @. (
    1
    + 0.17699 * (m14_x2_diff)
    - 0.50447 * (m14_x2_diff)^2
    - 0.02427 * (m14_x2_diff)^3
    + 0.72085 * (m14_x2_diff)^4
    + 0.01979 * (m14_x2_diff)^5
    - 0.77530 * (m14_x2_diff)^6
    + 0.32999 * (m14_x2_diff)^7
    + (0.0, 0.0, -0.011, 0.0, 0.0)
)
const m14_b2v = @. (
      1.41338 * (m14_x2 - 1.82)
    + 2.28305 * (m14_x2 - 1.82)^2
    + 1.07233 * (m14_x2 - 1.82)^3
    - 5.38434 * (m14_x2 - 1.82)^4
    - 0.62251 * (m14_x2 - 1.82)^5
    + 5.30260 * (m14_x2 - 1.82)^6
    - 2.09002 * (m14_x2 - 1.82)^7
    + (0.0, 0.0, 0.091, 0.0, 0.0)
)
const m14_a3v = @. (
      1.752
    - 0.316 * m14_x3
    - 0.104 / ((m14_x3 - 4.67)^2 + 0.341)
    + (0.442, 0.341, 0.130, 0.020, 0.000)
)
const m14_b3v = @. (
     -3.090
    + 1.825 * m14_x3
    + 1.206 / ((m14_x3 - 4.62)^2 + 0.263)
    - (1.256, 1.021, 0.416, 0.064, 0.000)
)
const m14_xn = collect((m14_x1..., m14_x2..., m14_x3...))
const m14_anv = collect((m14_a1v..., m14_a2v..., m14_a3v...))
const m14_bnv = collect((m14_b1v..., m14_b2v..., m14_b3v...))
const m14_a_spl = BSK.interpolate(m14_xn, m14_anv, BSK.BSplineOrder(4), BSK.Natural())
const m14_b_spl = BSK.interpolate(m14_xn, m14_bnv, BSK.BSplineOrder(4), BSK.Natural())

"""
    M14(;Rv=3.1)

Maiz Apellaniz et al (2014) Milky Way & LMC R(V) dependent model.

Returns A(λ)/A(V) at the given wavelength relative to the
extinction. The published UV extinction curve is identical to Clayton,
Cardelli, and Mathis (1989, CCM). Forcing the optical section to match smoothly
with CCM introduces a non-physical feature at high values of R5495 around 3.9
inverse microns; see section 5 in Maiz Apellaniz et al. (2014) for more
discussion. For that reason, we provide the M14 model only through 3.3 inverse
microns, the limit of the optical in CCM. Outside of that range this will
return 0. Rv is the selective extinction and is valid over [2, 6]. A typical
value for the Milky Way is 3.1. R5495 = A(5485)/E(4405-5495) Spectral
equivalent to photometric R(V).

# References
[Maiz Apellaniz et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014A%26A...564A..63M/)
"""
Parameters.@with_kw struct M14 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::M14)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return m14_invum(x, law.Rv)
end

bounds(::Type{M14}) = (3030.3030303030305, 33333.333333333336)

"""
    DustExtinction.m14_invum(x, Rv)

The algorithm used for the [`M14`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function m14_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 3.3)
        throw(DomainError(
            x, "out of bounds of M14, support is over $(bounds(M14)) angstrom"))
    end

    a = zero(x)
    b = zero(x)
    if x < m14_xi1
        # Infrared
        a = 0.574 * x^1.61
        b = -0.527 * x^1.61
    elseif x < m14_xi3
        # Optical
        a = m14_a_spl(x)
        b = m14_b_spl(x)
    end

    return a + b / Rv
end
