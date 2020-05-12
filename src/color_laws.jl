# Convenience function for wavelength conversion
@inline aa_to_invum(wave::Real) = 10000 / wave

# --------------------------------------------------------------------------------

# Optical coefficients
const ccm89_ca = [1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.7753, 0.32999, 0.0]
const ccm89_cb = [0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.3026, -2.09002, 0.0]

const od94_ca = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
const od94_cb = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]

"""
    CCM89(;Rv=3.1)

Clayton, Cardelli and Mathis (1989) dust law.

Returns E(B-V) in magnitudes at the given wavelength relative to the extinction
at 5494.5 Å. The default support is [1000, 33333]. Outside of that range this will return 0. `Rv` is the selective extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.

# References
[Clayton,Cardelli and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C)
"""
@with_kw struct CCM89 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::CCM89)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return ccm89_invum(x, law.Rv, ccm89_ca, ccm89_cb)
end

bounds(::Type{CCM89}) = (1000, 33333)


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

"""
    DustExtinction.ccm89_invum(x, Rv, c_a, c_b)

The algorithm used for the [`CCM89`](@ref) extinction law, given inverse microns, Rv, and a set of coefficients for use in the optical (only difference between ccm89 and od94). For more information, seek the original paper.
"""
function ccm89_invum(x::Real, Rv::Real, c_a::Vector{<:Real}, c_b::Vector{<:Real})
    if x < 0.3
        error("out of bounds of CCM89, support is over $(bounds(CCM89)) angstrom")
    elseif x < 1.1  # Near IR
        y = x^1.61
        a = 0.574y
        b = -0.527y
    elseif x < 3.3  # Optical
        y = x - 1.82
        yn = 1.0
        a = @evalpoly y c_a[1] c_a[2] c_a[3] c_a[4] c_a[5] c_a[6] c_a[7] c_a[8] c_a[9]
        b = @evalpoly y c_b[1] c_b[2] c_b[3] c_b[4] c_b[5] c_b[6] c_b[7] c_b[8] c_b[9]
    elseif x < 8.0  # NUV
        a =  1.752 - 0.316x - (0.104 / ((x - 4.67)^2 + 0.341))
        b = -3.090 + 1.825x + (1.206 / ((x - 4.62)^2 + 0.263))
        if x > 5.9 # Far NUV
            y = x - 5.9
            a += @evalpoly y 0.0 0.0 -0.04473 -0.009779
            b += @evalpoly y 0.0 0.0 0.213 0.1207
        end
    elseif x ≤ 10.0 # FUV
        y = x - 8.0
        a = @evalpoly y -1.073 -0.628 0.137 -0.07
        b = @evalpoly y 13.67 4.257 -0.42 0.374
    else
        error("out of bounds of CCM89, support is over $(bounds(CCM89)) angstrom")
    end
    return a + b / Rv
end

# --------------------------------------------------------------------------------

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
    return cal00_invum(x, law.Rv)
end

bounds(::Type{CAL00}) = (1200, 22000)

"""
    DustExtinction.cal00_invum(x, Rv)

The algorithm used for the [`CAL00`](@ref) extinction law, given inverse microns and Rv. For more information, seek the original paper.
"""
function cal00_invum(x::Real, Rv::Real)
    if x > 1 / 0.12
        error("out of bounds of CAL00, support is over $(bounds(CAL00)) angstrom")
    elseif x > 1 / 0.63
        k = @evalpoly x -2.156 1.509 -0.198 0.011
    elseif x > 1 / 2.2
        k = @evalpoly x -1.857 1.040
    else
        error("out of bounds of CAL00, support is over $(bounds(CAL00)) angstrom")
    end
    return 1.0 + 2.659 * k / Rv
end


"""
    GCC09(;Rv=3.1)

Gordon, Cartledge, & Clayton (2009) dust law.

This model applies to the UV spectral region all the way to 909.09 Å.
This model was not derived for the optical or NIR.
"""
@with_kw struct GCC09 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::GCC09)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return gcc09_invum(x, law.Rv)
end

bounds(::Type{GCC09}) = (909.09, 3030.3)

"""
    DustExtinction.gcc09_invum(x, Rv)

The algorithm used for the [`GCC09`](@ref) extinction law, given inverse microns and Rv. For more information, seek the original paper.
"""
function gcc09_invum(x::Real, Rv::Real)
    if 3.3 ≤ x ≤ 11.0  # NUV
        a = 1.894 - 0.373 * x - 0.0101 / ((x - 4.57)^2 + 0.0384)
        b = -3.490 + 2.057 * x + 0.706 / ((x - 4.59)^2 + 0.169)
    else # out of bounds
        error("out of bounds of GCC09, support is over $(bounds(GCC09)) angstrom")
    end
    if 5.9 ≤ x ≤ 11.0  # far-NUV
        y = x - 5.9
        a += @evalpoly y 0.0 0.0 -0.110 -0.0100
        b += @evalpoly y 0.0 0.0 0.531 0.0544
    end

    return a + b / Rv
end


"""
    VCG04(;Rv=3.1)

Valencic, Clayton, & Gordon (2004) dust law.

This model applies to the UV spectral region all the way to 912 Å.
This model was not derived for the optical or NIR.
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

The algorithm used for the [`VCG04`](@ref) extinction law, given inverse microns and Rv. For more information, seek the original paper.
"""
function vcg04_invum(x::Real, Rv::Real)
    if 3.3 ≤ x ≤ 8.0  # NUV
        a = 1.808 - 0.215 * x - 0.134 / ((x - 4.558)^2 + 0.566)
        b = -2.350 + 1.403 * x + 1.103 / ((x - 4.587)^2 + 0.263)
    else
        error("out of bounds of VCG04, support is over $(bounds(VCG04)) angstrom")
    end
    if 5.9 ≤ x ≤ 8.0  # far-NUV
        y = x - 5.9
        a += @evalpoly y 0.0 0.0 -0.0077 -0.0030
        b += @evalpoly y 0.0 0.0 0.2060 0.0550
    end

    return a + b / Rv
end
