# x value above which FM90 parametrization used
const f99_x_cutval_uv = aa_to_invum(2700)
# required UV points for spline interpolation
const f99_x_splineval_uv = aa_to_invum.((2700, 2600))

# Shape models used by F99 and F04
function _curve_F99_method(
        x,
        Rv,
        c1, c2, c3, c4,
        x0,
        gamma,
        optnir_axav_x, optnir_axav_y,
    )

    # add in required spline points, otherwise just spline points
    xuv = collect(f99_x_splineval_uv)
    if x >= f99_x_cutval_uv
        push!(xuv, x)
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
        spl = interpolate(collect(spline_x), collect(spline_y), BSplineOrder(4))
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
The default support is [1000, 33333] Å. Outside of that range this will
return 0. Rv is the selective extinction and is valid over [2, 6].
A typical value for the Milky Way is 3.1.

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

The algorithm used for the [`F99`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function f99_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 10.0)
        throw(DomainError(x, "out of bounds of F99, support is over $(bounds(F99)) angstrom"))
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

Returns E(B-V) in magnitudes at the given wavelength relative to the extinction.
This model applies to the UV and optical to NIR spectral range.
The default support is [1000, 33333] Å. Outside of that range this will return 0.
Rv is the selective extinction and is valid over [2, 6]. A typical value for
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

The algorithm used for the [`F04`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function f04_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 10.0)
        throw(DomainError(x, "out of bounds of F04, support is over $(bounds(F04)) angstrom"))
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
    nir_axebv_y = @. (0.63 * Rv - 0.84) * f04_nir_axav_x^1.84

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
The default support is [1149, 33333] Å. Outside of that range this will
return 0. Rv is the selective extinction and is valid over [2, 6].
A typical value for the Milky Way is 3.1.

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

bounds(::Type{F19}) = (1149.4252873563219, 33333.333333333336)

"""
    DustExtinction.f19_invum(x, Rv)

The algorithm used for the [`F19`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function f19_invum(x::Real, Rv::Real)
    # read and unpack tabulated data
    # using type annotations so `interpolate` can be inferred
    # TODO: Avoid using readdlm to avoid this issue altogether
    data_x::Vector{Float64}, data_k::Vector{Float64}, data_delta_k::Vector{Float64}, data_sigma_k::Vector{Float64} = let data = readdlm(joinpath(datadep"F19", "F19_tabulated.dat"), skipstart=1)
        (data[:, i] for i in 1:4)
    end

    if !(0.3 <= x <= 8.7)
        throw(DomainError(
            x, "out of bounds of F19, support is over $(bounds(F19)) angstrom"))
    end

    # compute E(lambda-55)/E(B-55) on the tabulated x points
    k_rV_tab_x = @. data_k + data_delta_k * (Rv - 3.10) * 0.990

    # setup spline interpolation
    spl = interpolate(collect(data_x), collect(k_rV_tab_x), BSplineOrder(4))

    # use spline interpolation to evaluate the curve for the input x values
    k_rV = spl(x)

    # convert to A(x)/A(55) from E(x-55)/E(44-55)
    a_rV = k_rV / Rv + 1.0

    return a_rV
end
