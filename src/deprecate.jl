# Here are all the deprecations

# deprecations from 0.4 on
@deprecate ccm89(x::AbstractArray, r_v::Real = 3.1) ccm89.(x, r_v)
@deprecate od94(x::AbstractArray, r_v::Real = 3.1) od94.(x, r_v)

# deprecations from 0.7 on
@deprecate redden(f::Real, λ::Real, Av::Real; Rv = 3.1, law = ccm89) = redden(CCM89, λ, f; Av = Av, Rv = 3.1)
@deprecate redden(f::Quantity, λ::Quantity, Av::Real; Rv = 3.1, law = ccm89) = redden(CCM89, λ, f; Av = Av, Rv = 3.1)

@deprecate deredden(f::Real, λ::Real, Av::Real; Rv = 3.1, law = ccm89) = deredden(CCM89, λ, f; Av = Av, Rv = 3.1)
@deprecate deredden(f::Quantity, λ::Quantity, Av::Real; Rv = 3.1, law = ccm89) = deredden(CCM89, λ, f; Av = Av, Rv = 3.1)

@deprecate ccm89(wave, Rv = 3.1) = CCM89(Rv = Rv)(wave)
@deprecate od94(wave, Rv = 3.1) = OD94(Rv = Rv)(wave)
@deprecate cal00(wave, Rv = 4.05) = CAL00(Rv = Rv)(wave)
