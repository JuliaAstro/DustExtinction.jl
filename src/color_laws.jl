using Dierckx

# Convenience function for wavelength conversion
@inline aa_to_invum(wave::Real) = 10000 / wave

# --------------------------------------------------------------------------------

# Optical coefficients
const ccm89_ca = [1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.7753, 0.32999, 0.0]
const ccm89_cb = [0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.3026, -2.09002, 0.0]

const od94_ca = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
const od94_cb = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]

# x value above which FM90 parametrization used
const f99_x_cutval_uv = aa_to_invum(2700)

# required UV points for spline interpolation
const f99_x_splineval_uv = aa_to_invum.((2700, 2600))

# spline points
const f99_optnir_axav_x = aa_to_invum.((26500, 12200, 6000, 5470, 4670, 4110))
const f99_opt_axebv_y_params = (
    (-0.426, +1.0044, +0.00000),
    (-0.050, +1.0016, +0.00000),
    (+0.701, +1.0016, +0.00000),
    (+1.208, +1.0032, -0.00033),
)
const f99_nir_axebv_y_params = @. (0.265, 0.829) / 3.1

# c1-c2 correlation terms
const f99_Rv_params = (-0.824, 4.717)
const f99_c1_params = (2.030, - 3.007)
const f99_c3 = 3.23
const f99_c4 = 0.41
const f99_x0 = 4.596
const f99_gamma = 0.99

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

"""
    GCC09(;Rv=3.1)

Gordon, Cartledge, & Clayton (2009) dust law.

This model applies to the UV spectral region all the way to 909.09 Å.
This model was not derived for the optical or NIR.

# References
[Gordon, Cartledge, & Clayton (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...705.1320G/)
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

# Shape models used by F99 and F04
function _curve_F99_method(
    x,
    Rv,
    c1,
    c2,
    c3,
    c4,
    x0,
    gamma,
    f99_optnir_axav_x,
    optnir_axav_y,
    )

    # add in required spline points, otherwise just spline points
    if x >= f99_x_cutval_uv
        xuv = (f99_x_splineval_uv..., x)
    else
        xuv = f99_x_splineval_uv
    end

    # FM90 model and values
    fm90_model = FM90(c1=c1, c2=c2, c3=c3, c4=c4, x0=x0, gamma=gamma)
    # evaluate model and get results in A(x)/A(V)
    axav_fm90 = @. fm90_model(aa_to_invum(xuv)) / Rv + 1

    # ignore the spline points
    if x >= f99_x_cutval_uv
        axav = last(axav_fm90)
    else
        # **Optical Portion**
        # using cubic spline anchored in UV, optical, and IR

        # optical/NIR points in input x

        # save spline points
        y_splineval_uv = axav_fm90

        # spline points
        x_splineval_optir = (0.0, f99_optnir_axav_x...)

        # determine optical/IR values at spline points
        y_splineval_optir = (0.0, optnir_axav_y...)
        spline_x = (x_splineval_optir..., f99_x_splineval_uv...)
        spline_y = (y_splineval_optir..., y_splineval_uv...)
        spl = Spline1D(collect(spline_x), collect(spline_y), k=3)
        axav = spl(x)
    end

    # return A(x)/A(V)
    return axav
end

"""
    F99(;Rv=3.1)

Fitzpatrick (1999) dust law.

Returns E(B-V) in magnitudes at the given wavelength relative to the
extinction. This model applies to the UV and optical to NIR spectral range.
The default support is [1000, 33333] Å. Outside of that range this will return
0. Rv is the selective extinction and is valid over [2, 6]. A typical value for
the Milky Way is 3.1.

# References
[Fitzpatrick (1999)](https://ui.adsabs.harvard.edu/abs/1999PASP..111...63F/)
"""
@with_kw struct F99 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::F99)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return f99_invum(x, law.Rv)
end

bounds(::Type{F99}) = (1000.0, 33333.3)

"""
    DustExtinction.f99_invum(x, Rv)

The algorithm used for the [`F99`](@ref) extinction law, given inverse microns and Rv. For more information, seek the original paper.
"""
function f99_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 10.0)
        error("out of bounds of F99, support is over $(bounds(F99)) angstrom")
    end

    # terms depending on Rv
    c2 = @evalpoly (1. / Rv) f99_Rv_params...
    # original F99 c1-c2 correlation
    c1 = @evalpoly c2 f99_c1_params...

    # determine optical/IR values at spline points
    #    Final optical spline point has a leading "-1.208" in Table 4
    #    of F99, but that does not reproduce Table 3.
    #    Additional indication that this is not correct is from
    #    fm_unred.pro
    #    which is based on FMRCURVE.pro distributed by Fitzpatrick.
    #    --> confirmation needed?
    #
    #    Also, fm_unred.pro has different coeff and # of terms,
    #    but later work does not include these terms
    #    --> check with Fitzpatrick?
    opt_axebv_y = evalpoly.(Rv, f99_opt_axebv_y_params)
    nir_axebv_y = @. f99_nir_axebv_y_params * Rv
    optnir_axebv_y = @. (nir_axebv_y..., opt_axebv_y...) / Rv

    return _curve_F99_method(
            x,
            Rv,
            c1,
            c2,
            f99_c3,
            f99_c4,
            f99_x0,
            f99_gamma,
            f99_optnir_axav_x,
            optnir_axebv_y,
        )
end
