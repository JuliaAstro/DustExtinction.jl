module DynamicQuantitiesExt

import DustExtinction: ExtinctionLaw, redden, deredden, SFD98Map, CSFDMap
using DynamicQuantities: UnionAbstractQuantity, ustrip, @u_str

# PARITY NOTE: this extension mirrors `UnitfulExt` for wavelength/angle inputs,
# but intentionally differs on the *output*. `UnitfulAstro` provides a
# logarithmic `mag` unit, so `UnitfulExt` returns `<value> * u"mag"`.
# `DynamicQuantities.jl` has no logarithmic-unit machinery, and forcing a
# stand-in dimension would misrepresent the quantity. We therefore return the
# extinction as a plain (dimensionless) number in magnitudes. Inputs are still
# converted to the expected units first, so the two extensions accept the same
# arguments and agree numerically.

# Evaluate a law at a wavelength `Quantity`, converting it to angstrom.
# Returns plain magnitudes (see PARITY NOTE above).
(l::ExtinctionLaw)(wavelength::UnionAbstractQuantity) = l(ustrip(u"Å", wavelength))

# Reddening / de-reddening with unitful wavelengths (and optionally unitful flux).
# When `flux` carries units the extinction factor is applied directly so the
# flux units are preserved (mirroring the Unitful behaviour of `flux * mag`).
redden(law::ExtinctionLaw, wave::UnionAbstractQuantity, flux::Real; Av = 1) =
    redden(law, ustrip(u"Å", wave), flux; Av)
redden(law::ExtinctionLaw, wave::UnionAbstractQuantity, flux::UnionAbstractQuantity; Av = 1) =
    flux * 10^(-0.4 * Av * law(ustrip(u"Å", wave)))
deredden(law::ExtinctionLaw, wave::UnionAbstractQuantity, flux::Real; Av = 1) =
    deredden(law, ustrip(u"Å", wave), flux; Av)
deredden(law::ExtinctionLaw, wave::UnionAbstractQuantity, flux::UnionAbstractQuantity; Av = 1) =
    flux / 10^(-0.4 * Av * law(ustrip(u"Å", wave)))

# Dust maps: angular coordinates given as quantities are converted to radians.
# Output is plain magnitudes (see PARITY NOTE above), whereas `UnitfulExt`
# returns `UnitfulAstro.mag`.
function (dustmap::SFD98Map)(l::UnionAbstractQuantity, b::UnionAbstractQuantity)
    return dustmap(ustrip(u"rad", l), ustrip(u"rad", b))
end

function (dustmap::CSFDMap)(l::UnionAbstractQuantity, b::UnionAbstractQuantity)
    return dustmap(ustrip(u"rad", l), ustrip(u"rad", b))
end

end # module
