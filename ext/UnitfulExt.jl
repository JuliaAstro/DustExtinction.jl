module UnitfulExt

import DustExtinction: ExtinctionLaw, redden, deredden, SFD98Map, CSFDMap, CCM89
import Unitful as U
using UnitfulAstro: UnitfulAstro

# Evaluate a law at a wavelength `Quantity`. The wavelength is converted to
# angstrom and the returned extinction is given as `UnitfulAstro.mag` (a
# logarithmic unit). The `DynamicQuantitiesExt` analog returns plain magnitudes
# instead, since DynamicQuantities.jl has no logarithmic-unit support.
(l::ExtinctionLaw)(wavelength::U.Quantity) = l(U.ustrip(U.u"Å", wavelength)) * U.u"mag"

# Reddening / de-reddening with unitful wavelengths (and optionally unitful flux)
redden(law::ExtinctionLaw, wave::U.Quantity, flux::Real; Av = 1) = redden(law, U.ustrip(U.u"Å", wave), flux; Av)
redden(law::ExtinctionLaw, wave::U.Quantity, flux::U.Quantity; Av = 1) = flux * (Av * law(wave))
deredden(law::ExtinctionLaw, wave::U.Quantity, flux::Real; Av = 1) = deredden(law, U.ustrip(U.u"Å", wave), flux; Av)
deredden(law::ExtinctionLaw, wave::U.Quantity, flux::U.Quantity; Av = 1) = flux / (Av * law(wave))

# Dust maps: galactic coordinates given as angular `Quantity`s are converted to
# radians and the returned extinction is given as `UnitfulAstro.mag`.
function (dustmap::SFD98Map)(l::U.Quantity, b::U.Quantity)
    l_ = U.ustrip(U.u"rad", l)
    b_ = U.ustrip(U.u"rad", b)
    return dustmap(l_, b_) * U.u"mag"
end

function (dustmap::CSFDMap)(l::U.Quantity, b::U.Quantity)
    l_ = U.ustrip(U.u"rad", l)
    b_ = U.ustrip(U.u"rad", b)
    return dustmap(l_, b_) * U.u"mag"
end

# Deprecations (moved here from src/deprecate.jl along with Unitful support)
@deprecate redden(f::U.Quantity, λ::U.Quantity, Av::Real; Rv = 3.1, law = CCM89) redden(law, λ, f; Av = Av, Rv = 3.1)
@deprecate deredden(f::U.Quantity, λ::U.Quantity, Av::Real; Rv = 3.1, law = CCM89) deredden(law, λ, f; Av = Av, Rv = 3.1)

end # module
