

"""
    FM90(;c1=0.10, c2=0.70, c3=3.23, c4=0.41, x0=4.60, gamma=0.9)
    FM90(coeffs, x0=4.60, gamma=0.9)

Fitzpatrick & Massa (1990) generative model for ultraviolet dust extinction. The default values are the published values for the Milky Way average.

## Parameters
* `c1` - y-intercept of linear term
* `c2` - slope of liner term
* `c3` - amplitude of 2175 Å bump
* `c4` - amplitude of FUV rise
* `x0` - centroid of 2175 Å bump
* `gamma` - width of 2175 Å bump

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the
returned value will be `UnitfulAstro.mag`.

## Examples
```jldoctest
julia> model = FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99);

julia> model(1500)
5.2521258452800135

julia> FM90()(1500)
5.152125845280013

julia> FM90(c1=0.2, c2=0.7, c3=3.23, c4=0.41, x0=4.6, gamma=0.99).([1000, 1200, 1800])
3-element Array{Float64,1}:
 12.562237969522851
  7.769215017329513
  4.890128210972148

```
# Extended Help

The model has form ``c_1 + c_2x + c_3D(x; \\gamma, x_0) + c_4 F(x)`` where ``x`` is the wavenumber in inverse microns, ``D(x)``
is a Drude profile (modified Lorentzian) used to model the 2175 Å bump with the scale-free parameters ``x_0`` (central wavenumber)
and ``\\gamma`` (damping coefficient), and ``F(x)``, a piecewise function for the far-UV. Note that the coefficients will change
the overall normalization, possibly changing the expected behavior of reddening via the parameter ``A_V``.

## References
[Fitzpatrick & Massa (1990)](https://ui.adsabs.harvard.edu/abs/1990ApJS...72..163F)
"""
@with_kw struct FM90{T<:Number}  <: ExtinctionLaw @deftype T
    c1 = 0.10
    c2 = 0.70
    c3 = 3.23
    c4 = 0.41
    x0 = 4.60
    gamma = 0.99
    @assert x0 ≥ 0 "`x0` must be ≥ 0, got $x0"
    @assert gamma ≥ 0 "`gamma` must be ≥ 0, got $gamma"
end

FM90(c1, c2, c3, c4, x0, gamma) = FM90(promote(c1, c2, c3, c4, x0, gamma)...)
FM90(coeffs, x0=4.60, gamma=0.99) = FM90(coeffs..., x0, gamma)

bounds(::Type{<:FM90}) = (912, 3200)

function (law::FM90)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))

    x = aa_to_invum(wave)
    exvebv = law.c1 + law.c2 * x

    x2 = x^2
    exvebv += law.c3 * (x2 / ((x2 - law.x0^2)^2 + x2 * (law.gamma^2)))

    if x >= 5.9
        y = x - 5.9
        exvebv += law.c4 * @evalpoly y 0 0 0.5392 0.05644
    end

    return exvebv
end


"""
    P92(BKG_amp=218.57142857142858, BKG_lambda=0.047, BKG_b=90.0, BKG_n=2.0, FUV_amp=18.545454545454547, FUV_lambda=0.07,
        FUV_b=4.0, FUV_n=6.5, NUV_amp=0.05961038961038961, NUV_lambda=0.22, NUV_b= -1.95, NUV_n=2.0, SIL1_amp=0.0026493506493506496,
        SIL1_lambda=9.7, SIL1_b=-1.95, SIL1_n=2.0, SIL2_amp=0.0026493506493506496, SIL2_lambda=18.0, SIL2_b= -1.8, SIL2_n=2.0,
        FIR_amp=0.015896103896103898, FIR_lambda=25.0, FIR_b=0.0, FIR_n=2.0)

Pei (1992) generative extinction model applicable from the extreme UV to far-IR.

## Parameters

Background Terms
* `BKG_amp` - amplitude
* `BKG_lambda` - central wavelength
* `BKG_b` - b coefficient
* `BKG_n` - n coefficient

Far-Ultraviolet Terms
* `FUV_amp` - amplitude
* `FUV_lambda` - central wavelength
* `FUV_b` - b coefficent
* `FUV_n` - n coefficient

Near-Ultraviolet (2175 Å) Terms
* `NUV_amp` - amplitude
* `NUV_lambda` - central wavelength
* `NUV_b` - b coefficent
* `NUV_n` - n coefficient

1st Silicate Feature (~10 micron) Terms
* `SIL1_amp` - amplitude
* `SIL1_lambda` - central wavelength
* `SIL1_b` - b coefficent
* `SIL1_n` - n coefficient

2nd Silicate Feature (~18 micron) Terms
* `SIL2_amp` - amplitude
* `SIL2_lambda` - central wavelength
* `SIL2_b` - b coefficient
* `SIL2_n` - n coefficient

Far-Infrared Terms
* `FIR_amp` - amplitude
* `FIR_lambda` - central wavelength
* `FIR_b` - b coefficent
* `FIR_n` - n coefficient

If `λ` is a `Unitful.Quantity` it will be automatically converted to Å and the returned value will be `UnitfulAstro.mag`.

## Examples
```jldoctest
julia> model = P92();

julia> model(1500)
2.396291891812002

julia> P92(FUV_b = 2.0).([1000, 2000, 3000])
3-element Array{Float64,1}:
 3.8390886792306187
 2.7304534614548697
 1.806181164464396

```

## Default Parameter Values

|Term |lambda|A|b|n|
|:---:|:---:|:---:|:---:|:---:|
|BKG  |0.047 |218.57142857142858   |90   |2  |
|FUV  |0.07  |18.545454545454547   |4.0  |6.5|
|NUV  |0.22  |0.05961038961038961  |-1.95|2.0|
|SIL1 |9.7   |0.0026493506493506496|-1.95|2.0|
|SIL2 |18.0  |0.0026493506493506496|-1.8 |2.0|
|FIR  |25.0  |0.015896103896103898 |0.0  |2.0|


## References
[Pei (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...395..130P)
"""
@with_kw struct P92{T<:Number} <: ExtinctionLaw @deftype T
    BKG_amp = 165.0 * (1 / 3.08 + 1)
    BKG_lambda = 0.047
    BKG_b = 90.0
    BKG_n = 2.0
    FUV_amp = 14.0 * (1 / 3.08 + 1)
    FUV_lambda = 0.07
    FUV_b = 4.0
    FUV_n = 6.5
    NUV_amp = 0.045 * (1 / 3.08 + 1)
    NUV_lambda = 0.22
    NUV_b = -1.95
    NUV_n = 2.0
    SIL1_amp = 0.002 * (1 / 3.08 + 1)
    SIL1_lambda = 9.7
    SIL1_b = -1.95
    SIL1_n = 2.0
    SIL2_amp = 0.002 * (1 / 3.08 + 1)
    SIL2_lambda = 18.0
    SIL2_b = -1.80
    SIL2_n = 2.0
    FIR_amp = 0.012 * (1 / 3.08 + 1)
    FIR_lambda = 25.0
    FIR_b = 0.00
    FIR_n = 2.0
    @assert BKG_amp ≥ 0 "`BKG_amp` must be ≥ 0, got $BKG_amp"
    @assert FUV_amp ≥ 0 "`FUV_amp` must be ≥ 0, got $FUV_amp"
    @assert 0.06 ≤ FUV_lambda ≤ 0.08 "`FUV_lambda` must be in between [0.06, 0.08], got $FUV_lambda"
    @assert NUV_amp ≥ 0 "`NUV_amp` must be ≥ 0, got $NUV_amp"
    @assert 0.20 ≤ NUV_lambda ≤ 0.24 "`NUV_lambda` must be in between [0.20, 0.24], got $NUV_lambda"
    @assert SIL1_amp ≥ 0 "`SIL1_amp` must be ≥ 0, got $SIL1_amp"
    @assert 7 ≤ SIL1_lambda ≤ 13 "`SIL1_lambda` must be in between [7, 13], got $SIL1_lambda"
    @assert SIL2_amp ≥ 0 "`SIL2_amp` must be ≥ 0, got $SIL2_amp"
    @assert 15 ≤ SIL2_lambda ≤ 21 "`SIL2_lambda` must be in between [15, 21], got $SIL2_lambda"
    @assert FIR_amp ≥ 0 "`FIR_amp` must be ≥ 0, got $FIR_amp"
    @assert 20 ≤ FIR_lambda ≤ 30 "`FIR_lambda` must be in between [20, 30], got $FIR_lambda"
end

P92(BKG_amp, BKG_lambda, BKG_b, BKG_n, FUV_amp, FUV_lambda, FUV_b, FUV_n,
    NUV_amp, NUV_lambda, NUV_b, NUV_n, SIL1_amp, SIL1_lambda, SIL1_b,
    SIL1_n, SIL2_amp, SIL2_lambda, SIL2_b, SIL2_n, FIR_amp, FIR_lambda,
    FIR_b, FIR_n) =
    P92(promote(BKG_amp, BKG_lambda, BKG_b, BKG_n, FUV_amp, FUV_lambda, FUV_b, FUV_n,
                NUV_amp, NUV_lambda, NUV_b, NUV_n, SIL1_amp, SIL1_lambda, SIL1_b,
                SIL1_n, SIL2_amp, SIL2_lambda, SIL2_b, SIL2_n, FIR_amp, FIR_lambda,
                FIR_b, FIR_n)...)

bounds(::Type{<:P92}) = (10, 10000000)

function (law::P92)(wave::T) where T
    checkbounds(law, wave) || return zero(float(T))

    x = aa_to_invum(wave)
    lam = 1.0 / x # wavelength is in microns

    axav = _p92_single_term(lam, law.BKG_amp, law.BKG_lambda, law.BKG_b, law.BKG_n)
    axav += _p92_single_term(lam, law.FUV_amp, law.FUV_lambda, law.FUV_b, law.FUV_n)
    axav += _p92_single_term(lam, law.NUV_amp, law.NUV_lambda, law.NUV_b, law.NUV_n)
    axav += _p92_single_term(lam, law.SIL1_amp, law.SIL1_lambda, law.SIL1_b, law.SIL1_n)
    axav += _p92_single_term(lam, law.SIL2_amp, law.SIL2_lambda, law.SIL2_b, law.SIL2_n)
    axav += _p92_single_term(lam, law.FIR_amp, law.FIR_lambda, law.FIR_b, law.FIR_n)

    return axav
end

# function for calculating a single P92 term
function _p92_single_term(x::Real, amplitude::Real, cen_wave::Real, b::Real, n::Real)
    l_norm = x / cen_wave
    return amplitude / (l_norm^n + inv(l_norm^n) + b)
end
