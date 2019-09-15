abstract type ColorLaw end

struct CCM89 <: ColorLaw
    RV::Real
end

function CCM89(RV=3.1)
    @assert 2 ≤ RV ≤ 6 "RV range spans [2.0, 6.0]"
    CCM89(RV)
end

struct OD94 <: ColorLaw
    RV::Real
end

function OD94(RV=3.1)
    @assert 2 ≤ RV ≤ 6 "RV range spans [2.0, 6.0]"
    OD94()
end

struct CAL00 <: ColorLaw
    RV::Real
end

function CAL00(RV=4.05)
    @assert 0 < RV "RV range spans [0, ∞)"
    CAL00(RV)
end 
