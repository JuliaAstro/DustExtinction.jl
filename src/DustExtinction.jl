module DustExtinction

using Unitful, UnitfulAstro

export extinct,
       extinct!,
       ccm89,
       cal00,
       od94,
       SFD98Map,
       ebv_galactic

# Convenience function for wavelength conversion
aa_to_invum(wave::Real) = 10000 / wave
aa_to_invum(wave::Quantity) = aa_to_invum(ustrip(u"angstrom", wave))


# Extinction Laws
include("ccm89.jl") # Also includes od94
include("cal00.jl")
include("SFD98Map.jl")

@deprecate ccm89(x::AbstractArray) ccm89.(x)
@deprecate od94(x::AbstractArray) od94.(x)

# Extinct function
"""
    extinct(f::Real, λ::Real, Av::Real; Rv::Real=3.1, law=ccm89)
    extinct(f::Quantity, λ::Quantity, Av::Real; Rv=3.1, law=ccm89)

Extinct the value `f` by the value calculated via the given law and total extinction value `Av`. By 
default we use `Rv=3.1` which is the Milky Way average selective attenuation.
"""
extinct(f::Real, λ::Real, Av::Real; Rv::Real = 3.1, law = ccm89) = f * 10^(-0.4 * Av * law(λ, Rv))
extinct(f::Quantity, λ::Quantity, Av::Real; Rv = 3.1, law = ccm89) = f * (Av * law(λ, Rv))

"""
    extinct!(f::AbstractArray{Real}, λ::AbstractArray{Real}, Av::Real; Rv::Real=3.1, law=ccm89)
    extinct!(f::AbstractArray{Quantity}, λ::AbstractArray{Quantity}, Av::Real; Rv=3.1, law=ccm89)

In-place version of [`extinct`](@ref) that only works with arrays.
"""
function extinct!(f::AbstractArray{<:Real}, λ::AbstractArray{<:Real}, Av::Real; Rv::Real = 3.1, law = ccm89) 
    @. f = f * 10^(-0.4 * Av * law(λ, Rv))
end
function extinct!(f::AbstractArray{<:Quantity}, λ::AbstractArray{<:Quantity}, Av::Real; Rv = 3.1, law = ccm89)
    @. f = f * (Av * law(λ, Rv))
end

end # module
