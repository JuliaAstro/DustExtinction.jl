module DustExtinction

# Convenience function for wavelength conversion
function aa_to_invum(wave::Real)
    return 1e4 / wave
end

# Extinction Laws
include("ccm89.jl") # Also includes od94
include("cal00.jl")
include("SFD98Map.jl")

end