using Polynomials

# Optical coefficients
const ccm89_ca = Poly([1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.7753, 0.32999])
const ccm89_cb = Poly([0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.3026, -2.09002])

const od94_ca = Poly([1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505])
const od94_cb = Poly([0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347])

function ccm89_invum(x::Real, r_v::Real, c_a::Poly{<:Real}, c_b::Poly{<:Real})
    if x < 0.3
        return 0.0x
    elseif x < 1.1  # Near IR
        y = x^1.61
        a = 0.574y
        b = -0.527y
    elseif x < 3.3  # Optical
        y = x - 1.82
        yn = 1.0
        a = c_a(y)
        b = c_b(y)
    elseif x < 8.0  # NUV
        a =  1.752 - 0.316x - (0.104 / ((x - 4.67)^2 + 0.341))
        b = -3.090 + 1.825x + (1.206 / ((x - 4.62)^2 + 0.263))
        if x > 5.9 # Far NUV
            y = x - 5.9
            a += Poly([0.0, 0.0, -0.04473, -0.009779])(y)
            b += Poly([0.0, 0.0, 0.213, 0.1207])(y)
        end
    elseif x ≤ 10.0 # FUV
        y = x - 8.0
        a = Poly([-1.073, -0.628, 0.137, -0.07])(y)
        b = Poly([13.67, 4.257, -0.42, 0.374])(y)
    else
        return 0.0x
    end

    return a + b / r_v
end

"""
    ccm89(wave::Real, r_v=3.1)
    ccm89(wave::Quantity, r_v=3.1)

Clayton, Cardelli and Mathis (1989) dust law. 

Returns the extinction in magnitudes at the given wavelength(s) `wave` (in Å) 
relative to the extinction at 5494.5 Å. Wavelength support is 1000. to 33333 Å. 
Outside of this range the returned value is 0. The parameter `r_v` changes the 
shape of the function.  A typical value for the Milky Way is 3.1.

# References
[[1]]
    (http://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C) Cardelli, Clayton 
    and Mathis (1989)
"""
function ccm89(wave::Real, r_v = 3.1)
    x = aa_to_invum(wave)
    return ccm89_invum(x, r_v, ccm89_ca, ccm89_cb)
end

ccm89(wave::Quantity, r_v::Real = 3.1) = ccm89(ustrip(u"angstrom", wave), r_v) * u"mag"

"""
    od94(wave::Real, r_v=3.1)
    od94(wave::Quantity, r_v=3.1)

O'Donnell (1994) dust law.

This is identical to the Clayton, Cardelli and Mathis (1989) dust law, except 
that different coefficients are used in the optical (3030.3 to 9090.9 Å). 
Returns the extinction in magnitudes at the given wavelength(s) `wave` (in Å) 
relative to the extinction at 5494.5 Å. Wavelength support is 1000 to 33333 Å. 
Outside of this range the returned value is 0. The parameter `r_v` changes 
the shape of the function. A typical value for the Milky Way is 3.1.

# References
[[1]]
    (http://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O) O'Donnell (1994)
"""
function od94(wave::Real, r_v = 3.1)
    x = aa_to_invum(wave)
    return ccm89_invum(x, r_v, od94_ca, od94_cb)
end
od94(wave::Quantity, r_v = 3.1) = ccm89(ustrip(u"angstrom", wave), r_v) * u"mag"
