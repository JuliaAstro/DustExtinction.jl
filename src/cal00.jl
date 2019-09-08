export cal00

function cal00_invum(x::Real, r_v::Real)

    if x > 1 / 0.12
        error("Wavelength out of range")
    elseif x > 1 / 0.63
        k = 2.659 * (((0.011x - 0.198)x + 1.509)x - 2.156)
    elseif x >  1 / 2.2
        k = 2.659 * (1.040x - 1.857)
    else
        error("Wavelength out of range")
    end

    return 1.0 + k / r_v

end

"""
    cal00(wave::Real, r_v::Real=3.1)
    cal00(wave::Quantity, r_v::Real=3.1)

Calzetti et al. (2000) Dust law.

Calculate the magnitudes for given wavelengths `wave` in Angstrom. Wavelength support is 0.12 to 2.2 microns (error will be thrown if out of this range). Accepts selective extinction `r_v` parameter with default set to Milky Way average of 3.1.

# References
[[1]](http://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C) Calzetti et al. (2000)

"""
function cal00(wave::Real, r_v::Real = 3.1)
    # Convert to inverse-um
    x = aa_to_invum.(wave)
    return cal00_invum(x, r_v)
end
cal00(wave::Quantity, r_v::Real = 3.1) = cal00(ustrip(u"angstrom", wave), r_v) * u"mag"
