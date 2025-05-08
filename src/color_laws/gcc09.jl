"""
    GCC09(;Rv=3.1)

Gordon, Cartledge, & Clayton (2009) dust law.

This model applies to the UV spectral region all the way to 909.09 Å.
This model was not derived for the optical or NIR.

# References
[Gordon, Cartledge, & Clayton (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...705.1320G/)
"""
Base.@kwdef struct GCC09 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::GCC09)(wave::T) where T <: Real
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return gcc09_invum(x, law.Rv)
end

bounds(::Type{GCC09}) = (909.0909090909091, 3030.3030303030305)

"""
    DustExtinction.gcc09_invum(x, Rv)

The algorithm used for the [`GCC09`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function gcc09_invum(x::Real, Rv::Real)
    if x < 3.3 || x > 11.0 # out of bounds
        throw(DomainError(x, "out of bounds of GCC09, support is over $(bounds(GCC09)) angstrom"))
    end
    # NUV
    a = 1.894 - 0.373 * x - 0.0101 / ((x - 4.57)^2 + 0.0384)
    b = -3.490 + 2.057 * x + 0.706 / ((x - 4.59)^2 + 0.169)
    if 5.9 ≤ x ≤ 11.0   # far-NUV
        y = x - 5.9
        a += @evalpoly y 0.0 0.0 -0.110 -0.0100
        b += @evalpoly y 0.0 0.0 0.531 0.0544
    end

    return a + b / Rv
end
