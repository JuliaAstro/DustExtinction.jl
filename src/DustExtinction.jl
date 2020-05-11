module DustExtinction

using Unitful
using UnitfulAstro
using DataDeps
using Parameters

export redden,
       deredden,
       CCM89,
       CAL00,
       OD94,
       GCC09,
       VCG04,
       SFD98Map,
       ebv_galactic


"""
    DustExtinction.ExtinctionLaw

The abstract super-type for dust extinction laws. See the extended help (`??DustExtinction.ExtinctionLaw` from the REPL) for more information about the interface.

## Extended Help
## Interface

Each extinction law implements the following methods
* [`DustExtinction.bounds(::ExtinctionLaw)`](@ref) - The bounds for the extinction law, as a `(min, max)` tuple in angstrom. If not implemented, will fallback to `(0, Inf)`
* [`(::ExtinctionLaw)(wavelength::Real)`](@ref) - the implmentation of the law, taking in angstrom and returning normalized extinction in astronomical magnitudes.

This is the bare-minimum required to use the law with [`redden`](@ref), [`deredden`](@ref), and the plotting recipes. Within the library we add support for [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) using code generation in `DustExtinction.jl/src/DustExtinction.jl`.
"""
abstract type ExtinctionLaw end

"""
    DustExtinction.bounds(::ExtinctionLaw)::Tuple

Get the natural wavelengths bounds for the extinction law, in angstrom
"""
bounds(::ExtinctionLaw) = (0, Inf)


"""
    redden(::ExtinctionLaw, wave, flux; Av=1)
    redden(::Type{ExtinctionLaw}, wave, flux; Av=1, kwargs...)

Redden the given `flux` using the given extinction law at the given wavelength.

If `wave` is `<:Real` then it is expected to be in angstrom and if it is `<:Unitful.Quantity` it will be automatically converted. `Av` is the total extinction value. The extinction law can be a constructed object or just a type. If it is just a type, `kwargs` will be passed to the constructor.

# Examples

```jldoctest
julia> wave = 3000; flux = 1000;

julia> redden(CCM89, wave, flux; Rv=3.1)

julia> redden(CCM89(Rv=3.1), wave, flux; Av=2)
35.11354215235764
```

# See Also
[`deredden`](@ref)
"""
redden(L::Type{<:ExtinctionLaw}, wave, flux; Av = 1, kwargs...) = redden(L(values(kwargs)...), wave, flux; Av = Av)
redden(law::ExtinctionLaw, wave::Real, flux; Av = 1) = flux * 10^(-0.4 * Av * law(wave))
redden(law::ExtinctionLaw, wave::Quantity, flux::Real; Av = 1) = redden(law, ustrip(u"Å", wave), flux; Av = Av)
redden(law::ExtinctionLaw, wave::Quantity, flux::Quantity; Av = 1) = flux * (Av * law(wave))

"""
    deredden(::ExtinctionLaw, wave, flux; Av=1)
    deredden(::Type{ExtinctionLaw}, wave, flux; Av=1, kwargs...)

Deredden the given `flux` using the given extinction law at the given wavelength.

If `wave` is `<:Real` then it is expected to be in angstrom and if it is `<:Unitful.Quantity` it will be automatically converted. `Av` is the total extinction value. The extinction law can be a constructed object or just a type. If it is just a type, `kwargs` will be passed to the constructor.

# Examples

```jldoctest
julia> wave = 3000:3005, flux = randn(size(wave));

julia> deredden(CCM89, wave, flux; Rv=3.1)

julia> deredden(CCM89(Rv=3.1), wave, flux; Av=2)
```

# See Also
[`redden`](@ref)
"""
deredden(L::Type{<:ExtinctionLaw}, wave, flux; Av = 1, kwargs...) = deredden(L(values(kwargs)...), wave, flux; Av = Av)
deredden(law::ExtinctionLaw, wave::Real, flux; Av = 1) = flux / 10^(-0.4 * Av * law(wave))
deredden(law::ExtinctionLaw, wave::Quantity, flux::Real; Av = 1) = deredden(law, ustrip(u"Å", wave), flux; Av = Av)
deredden(law::ExtinctionLaw, wave::Quantity, flux::Quantity; Av = 1) = flux / (Av * law(wave))

# --------------------------------------------------------------------------------
# bring in the support

include("deprecate.jl")
include("color_laws.jl")
include("dust_maps.jl")

# --------------------------------------------------------------------------------
# Here be codegen!

# generate unitful support for the following laws
# this can be removed when julia support is pinned to 1.3 or higher,
# at which point adding `(l::ExtinctionLaw)(wave)` is possible, until then
# using this code-gen does the trick but requires manually editing
# instead of providing support for all <: ExtinctionLaw
for law in [:CCM89, :OD94, :CAL00, :GCC09, :VCG04]
    @eval (l::$law)(wavelength::Quantity) = l(ustrip(u"Å", wavelength)) * u"mag"
end

# --------------------------------------------------------------------------------

function __init__()
    # register our data dependencies
    register(DataDep("sfd98_map",
    """
    SFD98 Galactic Dust Maps
    Website: https://sncosmo.github.io
    """,
    ["https://sncosmo.github.io/data/dust/SFD_dust_4096_ngp.fits",
    "https://sncosmo.github.io/data/dust/SFD_dust_4096_sgp.fits"],
    ["50b6aaad0b880762d0fd081177802dcc17c39d7044a410dd5649e2dfd0503e97",
    "84891a59054adab44a7be54051e4dcf0e66e3f13eee0d845ce3739242f553b83"]))
end

end # module
