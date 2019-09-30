using Polynomials, Unitful, UnitfulAstro

# Convenience function for wavelength conversion
aa_to_invum(wave::Real) = 10000 / wave
aa_to_invum(wave::Quantity) = aa_to_invum(ustrip(u"angstrom", wave))

#--------------------------------------------------------------------------------

# Optical coefficients
const ccm89_ca = Poly([1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.7753, 0.32999])
const ccm89_cb = Poly([0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.3026, -2.09002])

const od94_ca = Poly([1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505])
const od94_cb = Poly([0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347])

"""
    ccm89(λ::Real, Rv=3.1)
    ccm89(λ::Quantity, Rv=3.1)

Clayton, Cardelli and Mathis (1989) dust law. 

Returns E(B-V) in magnitudes at the given wavelength relative to the extinction 
at 5494.5 Å. `λ` is the wavelength in Å and has support over [1000, 33333]. 
Outside of that range this will return 0. `Rv` is the selective extinction 
and is valid over [2, 6]. A typical value for the Milky Way is 3.1

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the 
returned value will be `UnitfulAstro.mag`.

# References
[Clayton,Cardelli and Mathis (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C)
"""
function ccm89(λ::Real, Rv = 3.1)
    x = aa_to_invum(λ)
    return ccm89_invum(x, Rv, ccm89_ca, ccm89_cb)
end

ccm89(λ::Quantity, Rv::Real = 3.1) = ccm89(ustrip(u"angstrom", λ), Rv) * u"mag"

"""
    od94(λ::Real, Rv=3.1)
    od94(λ::Quantity, Rv=3.1)

O'Donnell (1994) dust law.

This is identical to the Clayton, Cardelli and Mathis (1989) dust law, except 
for different coefficients used in the optical (3030.3 Å to 9090.9 Å).

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and 
the returned value will be `UnitfulAstro.mag`.

# References
[O'Donnell (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O)

# See Also
[`ccm89`](@ref)
"""
function od94(λ::Real, Rv = 3.1)
    x = aa_to_invum(λ)
    return ccm89_invum(x, Rv, od94_ca, od94_cb)
end

od94(λ::Quantity, Rv = 3.1) = ccm89(ustrip(u"angstrom", λ), Rv) * u"mag"

function ccm89_invum(x::Real, Rv::Real, c_a::Poly{<:Real}, c_b::Poly{<:Real})
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

    return a + b / Rv
end

#--------------------------------------------------------------------------------

"""
    cal00(λ::Real, Rv=4.05)
    cal00(λ::Quantity, Rv=4.05)

Calzetti et al. (2000) Dust Law.

Returns E(B-V) in magnitudes at the given wavelength. `λ` is the wavelength in Å
 and has support over [1200, 22000]. Outside of that range this will return 0. 

Calzetti et al. (2000) developed a recipe for dereddening the spectra of 
galaxies where massive stars dominate the radiation output. They found the best
 fit value for such galaxies was 4.05±0.80.

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the 
returned value will be `UnitfulAstro.mag`.

# References
[Calzetti et al. (2000)](https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C)
"""
function cal00(λ::Real, Rv = 3.1)
    # Convert to inverse-um
    x = aa_to_invum.(λ)
    return cal00_invum(x, Rv)
end

cal00(λ::Quantity, Rv::Real = 3.1) = cal00(ustrip(u"angstrom", λ), Rv) * u"mag"

function cal00_invum(x::Real, Rv::Real)

    if x > 1 / 0.12
        return 0.0x
    elseif x > 1 / 0.63
        k = 2.659 * (((0.011x - 0.198)x + 1.509)x - 2.156)
    elseif x >  1 / 2.2
        k = 2.659 * (1.040x - 1.857)
    else
        return 0.0x
    end

    return 1.0 + k / Rv

end
