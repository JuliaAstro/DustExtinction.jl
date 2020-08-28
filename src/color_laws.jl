using Dierckx, DelimitedFiles

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
    if wave < 6300
        k = @evalpoly x -2.156 1.509 -0.198 0.011
    else
        k = @evalpoly x -1.857 1.040
    end
    return 1.0 + 2.659 * k / law.Rv
end

bounds(::Type{CAL00}) = (1200, 22000)

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

bounds(::Type{GCC09}) = (909.0909090909091, 3030.3030303030305)

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

# x value above which FM90 parametrization used
const f99_x_cutval_uv = aa_to_invum(2700)
# required UV points for spline interpolation
const f99_x_splineval_uv = aa_to_invum.((2700, 2600))

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
    optnir_axav_x,
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

        # save spline points
        y_splineval_uv = axav_fm90

        # spline points
        x_splineval_optir = (0.0, optnir_axav_x...)

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

# spline points
const f99_optnir_axav_x = aa_to_invum.((26500, 12200, 6000, 5470, 4670, 4110))
const f99_nir_axebv_y_params = @. (0.265, 0.829) / 3.1

# c1-c2 correlation terms
const f99_c3 = 3.23
const f99_c4 = 0.41
const f99_x0 = 4.596
const f99_gamma = 0.99

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
    c2 = @evalpoly (1. / Rv) -0.824 4.717
    # original F99 c1-c2 correlation
    c1 = @evalpoly c2 2.030 -3.007

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
    opt_axebv_y =
        (@evalpoly Rv -0.426 1.0044),
        (@evalpoly Rv -0.050 1.0016),
        (@evalpoly Rv  0.701 1.0016),
        (@evalpoly Rv  1.208 1.0032 -0.00033)

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

# spline points
const f04_opt_axav_x = aa_to_invum.((6000, 5470, 4670, 4110))
const f04_nir_axav_x = (0.50, 0.75, 1.0)

# **Use NIR spline x values in FM07, clipped to K band for now
const f04_optnir_axav_x = (f04_nir_axav_x..., f04_opt_axav_x...)

# c1-c2 correlation terms
const f04_c3 = 2.991
const f04_c4 = 0.319
const f04_x0 = 4.592
const f04_gamma = 0.922

"""
    F04(;Rv=3.1)
Fitzpatrick (2004) dust law.
Returns E(B-V) in magnitudes at the given wavelength relative to the
extinction. This model applies to the UV and optical to NIR spectral range.
The default support is [1000, 33333] Å. Outside of that range this will return
0. Rv is the selective extinction and is valid over [2, 6]. A typical value for
the Milky Way is 3.1.
Equivalent to the F99 model with an updated NIR Rv dependence
See also Fitzpatrick & Massa (2007, ApJ, 663, 320)
# References
[Fitzpatrick (2004)](https://ui.adsabs.harvard.edu/abs/2004ASPC..309...33F/)
"""
@with_kw struct F04 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::F04)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return f04_invum(x, law.Rv)
end

bounds(::Type{F04}) = (1000.0, 33333.3)

"""
    DustExtinction.f04_invum(x, Rv)
The algorithm used for the [`F04`](@ref) extinction law, given inverse microns and Rv. For more information, seek the original paper.
"""
function f04_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 10.0)
        error("out of bounds of F04, support is over $(bounds(F04)) angstrom")
    end

    # original F99 Rv dependence
    c2 = @evalpoly (1. / Rv) -0.824 4.717
    # updated F04 C1-C2 correlation
    c1 = @evalpoly c2 2.18 -2.91

    # **Keep optical spline points from F99:
    #    Final optical spline point has a leading "-1.208" in Table 4
    #    of F99, but that does not reproduce Table 3.
    #    Additional indication that this is not correct is from
    #    fm_unred.pro
    #    which is based on FMRCURVE.pro distributed by Fitzpatrick.
    #    --> confirmation needed?
    opt_axebv_y =
        (@evalpoly Rv -0.426 1.0044),
        (@evalpoly Rv -0.050 1.0016),
        (@evalpoly Rv  0.701 1.0016),
        (@evalpoly Rv  1.208 1.0032 -0.00033)

    # updated NIR curve from F04, note R dependence
    # Julia v1.0 workaround: https://github.com/JuliaAstro/DustExtinction.jl/pull/31#commitcomment-40605778
    nir_axebv_y_coeff = (0.63 * Rv - 0.84)
    nir_axebv_y = @. nir_axebv_y_coeff * f04_nir_axav_x^1.84

    optnir_axebv_y = @. (nir_axebv_y..., opt_axebv_y...) / Rv

    return _curve_F99_method(
            x,
            Rv,
            c1,
            c2,
            f04_c3,
            f04_c4,
            f04_x0,
            f04_gamma,
            f04_optnir_axav_x,
            optnir_axebv_y,
        )
end

"""
    F19(;Rv=3.1)

Fitzpatrick (2019) dust law.

Returns E(B-V) in magnitudes at the given wavelength relative to the
extinction. This model applies to the UV and optical to NIR spectral range.
The default support is [1149, 33333] Å. Outside of that range this will return
0. Rv is the selective extinction and is valid over [2, 6]. A typical value for
the Milky Way is 3.1.

Fitzpatrick, Massa, Gordon et al. (2019, ApJ, 886, 108) model. Based on a
sample of stars observed spectroscopically in the optical with HST/STIS.

# References
[Fitzpatrick (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...886..108F/)
"""
@with_kw struct F19 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::F19)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return f19_invum(x, law.Rv)
end

bounds(::Type{F19}) = (1149.4, 33333.3)

"""
    DustExtinction.f19_invum(x, Rv)

The algorithm used for the [`F19`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function f19_invum(x::Real, Rv::Real)
    # read and unpack tabulated data
    data_x, data_k, data_delta_k, data_sigma_k = let data = readdlm(joinpath(datadep"F19", "F19_tabulated.dat"), skipstart=1)
        (data[:, i] for i in 1:4)
    end

    if !(0.3 <= x <= 8.7)
        error("out of bounds of F19, support is over $(bounds(F19)) angstrom")
    end

    # compute E(lambda-55)/E(B-55) on the tabulated x points
    k_rV_tab_x = @. data_k + data_delta_k * (Rv - 3.10) * 0.990

    # setup spline interpolation
    spl = Spline1D(collect(data_x), collect(k_rV_tab_x), k=3)

    # use spline interpolation to evaluate the curve for the input x values
    k_rV = spl(x)

    # convert to A(x)/A(55) from E(x-55)/E(44-55)
    a_rV = k_rV / Rv + 1.0

    return a_rV
end

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
const m14_a_spl = Spline1D(m14_xn, m14_anv, bc="nearest")
const m14_b_spl = Spline1D(m14_xn, m14_bnv, bc="nearest")

"""
    M14(;Rv=3.1)

Maiz Apellaniz et al (2014) Milky Way & LMC R(V) dependent model.

Returns E(B-V) in magnitudes at the given wavelength relative to the
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
@with_kw struct M14 <: ExtinctionLaw
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

The algorithm used for the [`M14`](@ref) extinction law, given inverse microns and Rv. For more information, seek the original paper.
"""
function m14_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 3.3)
        error("out of bounds of M14, support is over $(bounds(M14)) angstrom")
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
