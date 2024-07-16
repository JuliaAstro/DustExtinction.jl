"""
    VCG04(;Rv=3.1)

Valencic, Clayton, & Gordon (2004) dust law.

This model applies to the UV spectral region all the way to 912 Å.
This model was not derived for the optical or NIR.

# References
[Valencic, Clayton, & Gordon (2004)](https://ui.adsabs.harvard.edu/abs/2004ApJ...616..912V/)
"""
@with_kw struct VCG04 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::VCG04)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return vcg04_invum(x, law.Rv)
end

bounds(::Type{VCG04}) = (1250, 3030.3)

"""
    DustExtinction.vcg04_invum(x, Rv)

The algorithm used for the [`VCG04`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function vcg04_invum(x::Real, Rv::Real)
    if x < 3.3 || x > 8.0
        throw(DomainError(x, "out of bounds of VCG04, support is over $(bounds(VCG04)) angstrom"))
    end

    # NUV
    a = 1.808 - 0.215 * x - 0.134 / ((x - 4.558)^2 + 0.566)
    b = -2.350 + 1.403 * x + 1.103 / ((x - 4.587)^2 + 0.263)
    if 5.9 ≤ x ≤ 8.0    # far-NUV
        y = x - 5.9
        a += @evalpoly y 0.0 0.0 -0.0077 -0.0030
        b += @evalpoly y 0.0 0.0 0.2060 0.0550
    end

    return a + b / Rv
end
