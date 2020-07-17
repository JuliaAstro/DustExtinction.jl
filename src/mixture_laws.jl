using DustExtinction: _curve_F99_method

const g03_obsdata_x = (
    0.455, 0.606, 0.800, 1.235, 1.538, 1.818, 2.273, 2.703,
    3.375, 3.625, 3.875, 4.125, 4.375, 4.625, 4.875, 5.125,
    5.375, 5.625, 5.875, 6.125, 6.375, 6.625, 6.875, 7.125,
    7.375, 7.625, 7.875, 8.125, 8.375, 8.625
)
const g03_obsdata_axav = (
    0.110, 0.169, 0.250, 0.567, 0.801, 1.000, 1.374, 1.672,
    2.000, 2.220, 2.428, 2.661, 2.947, 3.161, 3.293, 3.489,
    3.637, 3.866, 4.013, 4.243, 4.472, 4.776, 5.000, 5.272,
    5.575, 5.795, 6.074, 6.297, 6.436, 6.992
)

const g03_obsdata_tolerance = 6e-2

const g03_C1, g03_C2, g03_C3, g03_C4 = -4.959, 2.264, 0.389, 0.461
const g03_x0 = 4.6
const g03_gamma = 1.0

const g03_optnir_axav_x = 1.0 ./ (
    2.198, 1.65, 1.25, 0.81, 0.65, 0.55, 0.44, 0.37
)
# values at 2.198 and 1.25 changed to provide smooth interpolation
# as noted in Gordon et al. (2016, ApJ, 826, 104)
const g03_optnir_axav_y = (0.11, 0.169, 0.25, 0.567, 0.801, 1.00, 1.374, 1.672)

"""
    G03_SMCBar(;Rv=2.74) <Internal function>

Gordon et al. (2003) SMCBar Average Extinction Curve.

The observed A(lambda)/A(V) values at 2.198 and 1.25 microns were changed to
provide smooth interpolation as noted in Gordon et al. (2016, ApJ, 826, 104)

# Reference
[Gordon et al. (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...594..279G/)
"""
@with_kw struct G03_SMCBar <: DustExtinction.ExtinctionLaw
    Rv::Float64 = 2.74
    obsdata_x = g03_obsdata_x
    obsdata_axav = g03_obsdata_axav
    obsdata_tolerance::Float64 = g03_obsdata_tolerance
end

#G03_SMCBar(Rv) = G03_SMCBar(promote(Rv)...)
#G03_SMCBar(Rv=2.74) = G03_SMCBar(Rv)

bounds(::Type{<:G03_SMCBar}) = (1000.0, 33333.3)

function (law::G03_SMCBar)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return g03_invum(x, law.Rv)
end


function g03_invum(x::Real, Rv::Real)
    if !(0.3 <= x <= 10.0)
        error("out of bounds of G03, support is over $(bounds(G03)) angstrom")
    end

    # return A(x)/A(V)
    return _curve_F99_method(
        x,
        Rv,
        g03_C1,
        g03_C2,
        g03_C3,
        g03_C4,
        g03_x0,
        g03_gamma,
        g03_optnir_axav_x,
        g03_optnir_axav_y,
    )
end

"""
    G16(;Rv=3.1)

Gordon et al. (2016) Milky Way, LMC, & SMC R(V) and f_A dependent model

Returns E(B-V) in magnitudes at the given wavelength relative to the
extinction. This is mixture model between the F99 R(V) dependent model
(component A) and the [`G03_SMCBar`](@ref) model (component B) The default support is
[1000, 33333] Å. Outside of that range this will return 0. Rv is the selective
extinction and is valid over [2, 6]. A typical value for the Milky Way is 3.1.
By default, Rv=3.1 and f_A=1.0 is used.

# References
[Gordon et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...826..104G/)
"""
@with_kw struct G16{T<:Number} <: ExtinctionLaw @deftype T
    Rv::Float64 = 3.1
    f_A = 1.0
end

G16(Rv, f_A) = G16(promote(Rv, f_A)...)
#G16(Rv=3.1, f_A=3.1) = G16(Rv, f_A)

bounds(::Type{<:G16}) = (1000.0, 33333.3)

function (law::G16)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return g16_invum(x, law.Rv, law.f_A)
end

"""
    DustExtinction.g16_invum(x, Rv)
The algorithm used for the [`G16`](@ref) extinction law, given inverse microns
and Rv. For more information, seek the original paper.
"""
function g16_invum(x::Real, Rv::Real, f_A::Number)
    if !(0.3 <= x <= 10.0)
        error("out of bounds of G16, support is over $(bounds(G16)) angstrom")
    end

    # get the A component extinction model
    alav_A = f99_invum(x, Rv)

    # get the B component extinction model
    alav_B = g03_invum(x, Rv)

    # create the mixture model
    alav = @. f_A * alav_A + (1.0 - f_A) * alav_B

    return alav
end
