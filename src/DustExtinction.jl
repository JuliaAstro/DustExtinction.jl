module DustExtinction

using Unitful, UnitfulAstro

# Convenience function for wavelength conversion
aa_to_invum(wave::Real) = 10000 / wave
aa_to_invum(wave::Quantity) = aa_to_invum(ustrip(u"angstrom", wave))


# Extinction Laws
include("ccm89.jl") # Also includes od94
include("cal00.jl")
include("SFD98Map.jl")

@deprecate ccm89(x::AbstractArray) ccm89.(x::AbstractArray)
@deprecate od94(x::AbstractArray) od94.(x::AbstractArray)

end
