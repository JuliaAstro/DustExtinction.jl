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

# read and unpack tabulated data
# from https://raw.githubusercontent.com/karllark/dust_extinction/refs/heads/master/dust_extinction/data/F19_tabulated.dat
# TODO: Consider other places to store and track this dataset. See discussion in https://github.com/JuliaAstro/DustExtinction.jl/pull/54
const f99_x, f99_k, f99_delta_k, _ = eachcol([
 0.000  -3.020  -1.000   0.000
 0.455  -2.747  -0.842   0.000
 0.606  -2.528  -0.728   0.068
 0.800  -2.222  -0.531   0.070
 1.000  -1.757  -0.360   0.184
 1.100  -1.567  -0.284   0.100
 1.200  -1.300  -0.223   0.085
 1.250  -1.216  -0.198   0.072
 1.300  -1.070  -0.173   0.073
 1.350  -0.973  -0.150   0.061
 1.400  -0.868  -0.130   0.053
 1.450  -0.750  -0.110   0.052
 1.500  -0.629  -0.096   0.048
 1.550  -0.509  -0.081   0.041
 1.600  -0.407  -0.063   0.044
 1.650  -0.320  -0.048   0.037
 1.700  -0.221  -0.032   0.033
 1.750  -0.133  -0.017   0.030
 1.800  -0.048  -0.005   0.022
 1.818   0.000   0.000   0.022
 1.850   0.071   0.007   0.018
 1.900   0.188   0.013   0.021
 1.950   0.319   0.012   0.026
 2.000   0.438   0.010   0.032
 2.050   0.575   0.004   0.034
 2.100   0.665   0.003   0.027
 2.150   0.744   0.000   0.024
 2.200   0.838   0.002   0.017
 2.250   0.951   0.001   0.017
 2.273   1.000   0.000   0.016
 2.300   1.044  -0.000   0.020
 2.350   1.113   0.001   0.022
 2.400   1.181   0.001   0.027
 2.450   1.269  -0.002   0.034
 2.500   1.346   0.000   0.041
 2.550   1.405  -0.002   0.048
 2.600   1.476  -0.002   0.063
 2.650   1.558  -0.006   0.075
 2.700   1.632  -0.009   0.085
 2.750   1.723  -0.011   0.094
 2.800   1.791  -0.017   0.103
 2.850   1.869  -0.025   0.112
 2.900   1.948  -0.029   0.118
 2.950   2.009  -0.037   0.122
 3.000   2.090  -0.043   0.127
 3.100   2.253  -0.064   0.132
 3.200   2.408  -0.092   0.139
 3.300   2.565  -0.122   0.152
 3.400   2.746  -0.161   0.164
 3.500   2.933  -0.201   0.176
 3.600   3.124  -0.249   0.189
 3.700   3.328  -0.303   0.205
 3.800   3.550  -0.366   0.224
 3.900   3.815  -0.437   0.247
 4.000   4.139  -0.517   0.277
 4.100   4.534  -0.603   0.318
 4.200   5.012  -0.692   0.374
 4.300   5.560  -0.774   0.454
 4.400   6.118  -0.843   0.559
 4.500   6.565  -0.888   0.665
 4.600   6.767  -0.908   0.714
 4.700   6.681  -0.903   0.675
 4.800   6.394  -0.880   0.596
 4.900   6.038  -0.849   0.524
 5.000   5.704  -0.816   0.475
 5.100   5.432  -0.785   0.445
 5.200   5.226  -0.760   0.431
 5.300   5.078  -0.741   0.426
 5.400   4.978  -0.729   0.428
 5.500   4.913  -0.722   0.436
 5.600   4.877  -0.722   0.447
 5.700   4.862  -0.726   0.461
 5.800   4.864  -0.734   0.477
 5.900   4.879  -0.745   0.496
 6.000   4.904  -0.760   0.516
 6.100   4.938  -0.778   0.538
 6.200   4.982  -0.798   0.562
 6.300   5.038  -0.820   0.587
 6.400   5.105  -0.845   0.613
 6.500   5.181  -0.870   0.640
 6.600   5.266  -0.898   0.668
 6.700   5.359  -0.926   0.698
 6.800   5.460  -0.956   0.728
 6.900   5.569  -0.988   0.760
 7.000   5.684  -1.020   0.793
 7.100   5.805  -1.053   0.827
 7.200   5.933  -1.087   0.861
 7.300   6.067  -1.122   0.897
 7.400   6.207  -1.158   0.934
 7.500   6.352  -1.195   0.972
 7.600   6.502  -1.232   1.011
 7.700   6.657  -1.270   1.051
 7.800   6.817  -1.309   1.091
 7.900   6.981  -1.349   1.133
 8.000   7.150  -1.389   1.176
 8.100   7.323  -1.429   1.220
 8.200   7.500  -1.471   1.264
 8.300   7.681  -1.513   1.310
 8.400   7.866  -1.555   1.357
 8.500   8.054  -1.598   1.404
 8.600   8.246  -1.641   1.453
 8.700   8.441  -1.685   1.502
])

"""
    DustExtinction.f19_invum(x, Rv)

The algorithm used for the [`F19`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function f19_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 8.7)
        throw(DomainError(
            x, "out of bounds of F19, support is over $(bounds(F19)) angstrom"))
    end

    # compute E(lambda-55)/E(B-55) on the tabulated x points
    k_rV_tab_x = @. f99_k + f99_delta_k * (Rv - 3.10) * 0.990

    # setup spline interpolation
    spl = interpolate(f99_x, k_rV_tab_x, BSplineOrder(4))

    # use spline interpolation to evaluate the curve for the input x values
    k_rV = spl(x)

    # convert to A(x)/A(55) from E(x-55)/E(44-55)
    a_rV = k_rV / Rv + 1.0

    return a_rV
end
