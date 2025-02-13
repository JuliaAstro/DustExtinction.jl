using BSplineKit: interpolate, BSplineOrder, Natural
using DelimitedFiles

# Convenience function for wavelength conversion
@inline aa_to_invum(wave::Real) = 10_000 / wave

# --------------------------------------------------------------------------------


"""
    OD94(;Rv=3.1)

O'Donnell (1994) dust law.

This is identical to the Clayton, Cardelli and Mathis (1989) dust law, except
for different coefficients used in the optical (3030.3 Å to 9090.9 Å).

# References
[O'Donnell (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O)

# See Also
[`CCM89`](@ref)
"""
@with_kw struct OD94 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::OD94)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return ccm89_invum(x, law.Rv, od94_ca, od94_cb)
end

bounds(::Type{OD94}) = (1000, 33333)

# ------------------------------------------------------------------------------

"""
    CAL00(;Rv=4.05)

Calzetti et al. (2000) Dust Law.

Returns E(B-V) in magnitudes at the given wavelength. `λ` is the wavelength in Å
and has support over [1200, 22000]. Outside of that range this will return 0.

Calzetti et al. (2000) developed a recipe for dereddening the spectra of
galaxies where massive stars dominate the radiation output. They found the best
fit value for such galaxies was 4.05±0.80.

# References
[Calzetti et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C)
"""
@with_kw struct CAL00 <: ExtinctionLaw
    Rv::Float64 = 4.05
end
function (law::CAL00)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    if wave < 6300
        k = @evalpoly x -2.156 1.509 -0.198 0.011
    else
        k = @evalpoly x -1.857 1.040
    end
    return 1.0 + 2.659 * k / law.Rv
end

bounds(::Type{CAL00}) = (1200, 22000)

include("color_laws/ccm89.jl")
include("color_laws/fitzpatrick.jl")
include("color_laws/gcc09.jl")
include("color_laws/m14.jl")
include("color_laws/vcg04.jl")
