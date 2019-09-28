module DustExtinction

using Unitful, UnitfulAstro

export redden,
       deredden,
       ccm89,
       cal00,
       od94,
       SFD98Map,
       ebv_galactic


include("color_laws.jl")
include("dust_maps.jl")

@deprecate ccm89(x::AbstractArray, r_v::Real=3.1) ccm89.(x, r_v)
@deprecate od94(x::AbstractArray, r_v::Real=3.1) od94.(x, r_v)

# reddening functions
"""
    redden(f::Real, λ::Real, Av; RV=3.1, law=ccm89)
    redden(f::Quantity, λ::Quantity, Av; RV=3.1, law=ccm89)

Redden the value `f` by the value calculated via the given law and total 
extinction value `Av`. By default we use `RV=3.1` which is the Milky Way 
average selective attenuation. Note that λ should be in Angstrom if it is not 
a `Quantity`.
"""
redden(f::Real, λ::Real, Av::Real; RV = 3.1, law = ccm89) = f * 10^(-0.4 * Av * law(λ, RV))
redden(f::Quantity, λ::Quantity, Av::Real; RV = 3.1, law = ccm89) = f * (Av * law(λ, RV))

"""
    deredden(f::Real, λ::Real, Av; RV=3.1, law=ccm89)
    deredden(f::Quantity, λ::Quantity, Av; RV=3.1, law=ccm89)

Deredden the value `f` by the value calculated via the given law and total 
extinction value `Av`. By default we use `RV=3.1` which is the Milky Way 
average selective attenuation. Note that λ should be in Angstrom if it is not 
a `Quantity`.
"""
deredden(f::Real, λ::Real, Av::Real; RV = 3.1, law = ccm89) = f / 10^(-0.4 * Av * law(λ, RV))
deredden(f::Quantity, λ::Quantity, Av::Real; RV = 3.1, law = ccm89) = f / (Av * law(λ, RV))

end # module
