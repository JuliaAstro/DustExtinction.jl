module DustExtinction

using Unitful
using UnitfulAstro
using DataDeps
using Parameters

export redden,
       deredden,
       # Rv laws
       CCM89,
       CAL00,
       OD94,
       GCC09,
       VCG04,
       F99,
       F04,
       F19,
       M14,
       # Fittable laws
       FM90,
       # Mixture laws
       G16,
       G03_SMCBar,
       G03_LMCAve,
       # Dust maps
       SFD98Map,
       ebv_galactic,
       # Makie.jl recipe extension
       lplot,
       lplot!,
       dplot,
       dplot!

# Plot recipe stubs
function lplot end
function lplot! end
function dplot end
function dplot! end

"""
    abstract type DustExtinction.ExtinctionLaw

The abstract supertype for dust extinction laws. See the extended help
(`??DustExtinction.ExtinctionLaw` from the REPL) for more information about the interface.

# Extended Help

## Interface

Here's how to make a new extinction law, called `MyLaw`
* Create your struct. We strongly recommend using `Parameters.jl` to facilitate
  creating keyword argument constructors if your model is parameterized, which
  allows convenient usage with [`redden`](@ref) and [`deredden`](@ref).
  ```julia
  struct MyLaw <: DustExtinction.ExtinctionLaw end
  ```
* (Optional) Define the limits. This will default to `(0, Inf)`. Currently, this
  is used within the [`DustExtinction.checkbounds`](@ref) function and in the
  future will be used for plotting recipes.
  ```julia
  DustExtinction.bounds(::Type{<:MyLaw}) = (min, max)
  ```
* Define the law. You only need to provide one function which takes wavelength
  as angstrom. If your law is naturally written for inverse-micron, there is a
  helper function `aa_to_invum`.
  ```julia
  (::MyLaw)(wavelength::Real)
  ```
* (Optional) enable `Unitful.jl` support by adding this function. If you are
  building a new law within `DustExtinction.jl` you can add your law to the
  code-gen list inside `DustExtinction.jl/src/DustExtinction.jl`.
  ```julia
  (l::MyLaw)(wavelength::Unitful.Quantity) = l(ustrip(u"angstrom", wavelength)) * u"mag"
  ```
"""
abstract type ExtinctionLaw end

Base.broadcastable(law::ExtinctionLaw) = (law,)

"""
    DustExtinction.bounds(::ExtinctionLaw)::Tuple
    DustExtinction.bounds(::Type{<:ExtinctionLaw})::Tuple

Get the natural wavelengths bounds for the extinction law, in angstrom
"""
bounds(::E) where {E <: ExtinctionLaw} = bounds(E)
bounds(::Type{<:ExtinctionLaw}) = (0, Inf)

"""
    DustExtinction.checkbounds(::ExtinctionLaw, wavelength)::Bool
    DustExtinction.checkbounds(::Type{<:ExtinctionLaw, wavelength}::Bool

Helper function that uses [`DustExtinction.bounds`](@ref) to return whether the
given wavelength is in the support for the law.
"""
checkbounds(::E, wave) where {E <: ExtinctionLaw} = checkbounds(E, wave)
function checkbounds(E::Type{<:ExtinctionLaw}, wave)
    b = bounds(E)
    return b[1] ≤ wave ≤ b[2]
end


"""
    redden(::ExtinctionLaw, wave, flux; Av=1)
    redden(::Type{ExtinctionLaw}, wave, flux; Av=1, law_kwargs...)

Redden the given `flux` using the given extinction law at the given wavelength.

If `wave` is `<:Real` then it is expected to be in angstrom and if it is
`<:Unitful.Quantity` it will be automatically converted. `Av` is the total
extinction value. The extinction law can be a constructed struct or a `Type`.
If it is a `Type`, `law_kwargs` will be passed to the constructor.

# Examples

```jldoctest
julia> wave = 3000; flux = 1000;

julia> redden(CCM89, wave, flux; Rv=3.1)
187.38607779757183

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
    deredden(::Type{ExtinctionLaw}, wave, flux; Av=1, law_kwargs...)

Deredden the given `flux` using the given extinction law at the given wavelength.

If `wave` is `<:Real` then it is expected to be in angstrom and if it is
`<:Unitful.Quantity` it will be automatically converted. `Av` is the total
extinction value. The extinction law can be a constructed struct or a `Type`.
If it is a `Type`, `law_kwargs` will be passed to the constructor.

# Examples

```jldoctest
julia> wave = 3000; flux = 187.386;

julia> deredden(CCM89, wave, flux; Rv=3.1)
999.9995848273642

julia> deredden(CCM89(Rv=3.1), wave, flux; Av=2)
5336.573541539394
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
include("fittable_laws.jl")
include("mixture_laws.jl")

# --------------------------------------------------------------------------------
# Here be codegen!

# generate unitful support for the following laws
# this can be removed when julia support is pinned to 1.3 or higher,
# at which point adding `(l::ExtinctionLaw)(wave::Quantity)` is possible, until then
# using this code-gen does the trick but requires manually editing
# instead of providing support for all <: ExtinctionLaw
for law in [CCM89, OD94, CAL00, GCC09, VCG04, FM90, G16, G03_SMCBar, G03_LMCAve, F99, F04, F19, M14]
    (l::law)(wavelength::Quantity) = l(ustrip(u"Å", wavelength)) * u"mag"
end

# --------------------------------------------------------------------------------

function __init__()
    # register our data dependencies
    register(
        DataDep(
            "sfd98_map",
            """
            SFD98 Galactic Dust Maps
            Website: https://sncosmo.github.io
            """,
            ["https://sncosmo.github.io/data/dust/SFD_dust_4096_ngp.fits",
             "https://sncosmo.github.io/data/dust/SFD_dust_4096_sgp.fits"],
            ["50b6aaad0b880762d0fd081177802dcc17c39d7044a410dd5649e2dfd0503e97",
             "84891a59054adab44a7be54051e4dcf0e66e3f13eee0d845ce3739242f553b83"]
        )
    )
end

end # module
