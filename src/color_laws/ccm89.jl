# Optical coefficients
const ccm89_ca = [1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.7753, 0.32999, 0.0]
const ccm89_cb = [0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.3026, -2.09002, 0.0]

const od94_ca = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
const od94_cb = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]

"""
    CCM89(;Rv=3.1)

Clayton, Cardelli and Mathis (1989) dust law.

Returns A(λ)/A(V) at the given wavelength relative to the extinction
at 5494.5 Å. The default support is [1000, 33333]. Outside of that range this
will return 0. `Rv` is the selective extinction and is valid over [2, 6].
A typical value for the Milky Way is 3.1.

# References
[Clayton,Cardelli and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C)
"""
Base.@kwdef struct CCM89 <: ExtinctionLaw
    Rv::Float64 = 3.1
end

function (law::CCM89)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))
    x = aa_to_invum(wave)
    return ccm89_invum(x, law.Rv, ccm89_ca, ccm89_cb)
end

bounds(::Type{CCM89}) = (1000, 33333)

"""
    DustExtinction.ccm89_invum(x, Rv, c_a, c_b)

The algorithm used for the [`CCM89`](@ref) extinction law, given inverse microns,
Rv, and a set of coefficients for use in the optical (only difference between
ccm89 and od94). For more information, seek the original paper.
"""
function ccm89_invum(x::Real, Rv::Real, c_a::Vector{<:Real}, c_b::Vector{<:Real})
    if x < 0.3 || x > 10.0
        throw(DomainError(x, "out of bounds of CCM89, support is over $(bounds(CCM89)) angstrom"))
    end
    if x < 1.1      # Near IR
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
    else            # FUV
        y = x - 8.0
        a = @evalpoly y -1.073 -0.628 0.137 -0.07
        b = @evalpoly y 13.67 4.257 -0.42 0.374
    end
    return a + b / Rv
end
