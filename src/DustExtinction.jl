module DustExtinction

using Unitful, UnitfulAstro, DataDeps

export redden,
       deredden,
       ccm89,
       cal00,
       od94,
       gcc09,
       vcg04,
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

Redden the given `flux` by the extinction law at the given wavelength.

If `wave` is `<:Real` then it is expected to be in angstrom and if it is `<:Unitful.Quantity` it will be automatically converted. `Av` is the total extinction value. The extinction law can be a constructed object or just a type. If it is just a type, `kwargs` will be passed to the constructor.

# Examples

```jldoctest
julia> wave = 3000:3005, flux = randn(size(wave));

julia> redden(CCM89, wave, flux; Rv=3.1)

julia> redden(CCM89(Rv=3.1), wave, flux; Av=2)
```

# See Also
[`deredden`](@ref)
"""
redden(L::Type{<:ExtinctionLaw}, wave, flux; Av = 1, kwargs...) = redden(L(kwargs...), wave, flux; Av = Av)
redden(law::ExtinctionLaw, wave::Real, flux; Av = 1) = flux * 10^(-0.4 * Av * law(wave))
redden(law::ExtinctionLaw, wave::Quantity, flux; Av = 1) = flux * (Av * law(wave))

"""
    deredden(f::Real, λ::Real, Av; Rv=3.1, law=ccm89)
    deredden(f::Quantity, λ::Quantity, Av; Rv=3.1, law=ccm89)

Deredden the value `f` by the value calculated via the given law and total
extinction value `Av`. By default we use `Rv=3.1` which is the Milky Way
average selective attenuation. Note that λ should be in Angstrom if it is not
a `Quantity`.
"""
deredden(f::Real, λ::Real, Av::Real; Rv = 3.1, law = ccm89) = f / 10^(-0.4 * Av * law(λ, Rv))
deredden(f::Quantity, λ::Quantity, Av::Real; Rv = 3.1, law = ccm89) = f / (Av * law(λ, Rv))

# --------------------------------------------------------------------------------
# bring in the laws

include("deprecate.jl")
include("color_laws.jl")
include("dust_maps.jl")



# --------------------------------------------------------------------------------
# Here be codegen!

# generature unitful support for the following laws
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
