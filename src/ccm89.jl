export ccm89, od94

# Optical coefficients
const ccm89_ca = [1., 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530,
                  0.32999]
const ccm89_cb = [0., 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260,
                  -2.09002]
const od94_ca = [1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647,
                 -0.505]
const od94_cb = [0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805,
                 3.347]

function ccm89_invum(x::Real, r_v::Real, c_a::Vector{<:Real}, c_b::Vector{<:Real})
    a = 0.
    b = 0.
    if x < 0.3
        error("wavelength out of range")
    elseif x < 1.1  # Near IR
        y = x^1.61
        a = 0.574y
        b = -0.527y
    elseif x < 3.3  # Optical
        y = x - 1.82
        yn = 1.
        a = c_a[1]
        b = c_b[1]
        for i = 2:length(c_a)
            yn *= y
            a += c_a[i] * yn
            b += c_b[i] * yn
        end
    elseif x < 8.  # UV
        a =  1.752 - 0.316x - (0.104 / ((x - 4.67)^2 + 0.341))
        b = -3.090 + 1.825x + (1.206 / ((x - 4.62)^2 + 0.263))
        if x > 5.9
            y = x - 5.9
            y2 = y * y
            y3 = y2 * y
            a += -0.04473y2 - 0.009779y3
            b += 0.2130y2 + 0.1207y3
        end
    elseif x < 10.
        y = x - 8.
        y2 = y * y
        y3 = y2 * y
        a = -0.070y3 + 0.137y2 - 0.628y - 1.073
        b = 0.374y3 - 0.420y2 + 4.257y + 13.670
    else
        error("wavelength out of range")
    end

    a + b / r_v
end

"""
    ccm89(wave::Real, r_v::Real=3.1)
    ccm89(wave::Quantity, r_v::Real=3.1)

Clayton, Cardelli and Mathis (1989) dust law. 

Returns the extinction in magnitudes at the given wavelength(s) `wave` (in Angstroms),
relative to the extinction at 5494.5 Angstroms. The parameter `r_v`
changes the shape of the function.  A typical value for the Milky Way
is 3.1. An error is raised for wavelength values outside the range of
support, 1000. to 33333.33 Angstroms.

# References
[[1]]
    (http://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C) Cardelli, Clayton and Mathis (1989)
"""
function ccm89(wave::Real, r_v::Real = 3.1)
    x = aa_to_invum(wave)
    return ccm89_invum(x, r_v, ccm89_ca, ccm89_cb)
end

ccm89(wave::Quantity, r_v::Real = 3.1) = ccm89(ustrip(u"angstrom", wave), r_v) * u"mag"

"""
    od94(wave::Real, r_v::Real=3.1)
    od94(wave::Quantity, r_v::Real=3.1)

O'Donnell (1994) dust law.

This is identical to the Clayton, Cardelli
and Mathis (1989) dust law, except that different coefficients are
used in the optical (3030.3 to 9090.9 Angstroms). Returns the
extinction in magnitudes at the given wavelength(s) `wave` (in
Angstroms), relative to the extinction at 5494.5 Angstroms. The
parameter `r_v` changes the shape of the function.  A typical value
for the Milky Way is 3.1.  An error is raised for wavelength values
outside the range of support, 1000. to 33333.33 Angstroms.

# References
[[1]](http://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O) O'Donnell (1994)
"""
function od94(wave::Real, r_v::Real = 3.1)
    x = aa_to_invum(wave)
    return ccm89_invum(x, r_v, od94_ca, od94_cb)
end
od94(wave::Quantity, r_v::Real = 3.1) = ccm89(ustrip(u"angstrom", wave), r_v) * u"mag"
