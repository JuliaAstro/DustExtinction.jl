using FITSIO

import Base: show

# SFD98 Dust Maps

mutable struct SFD98Map
    mapdir::String
    ngp::ImageHDU
    ngp_size::Tuple{Int,Int}
    ngp_crpix1::Float64
    ngp_crpix2::Float64
    ngp_lam_scal::Float64
    sgp::ImageHDU
    sgp_size::Tuple{Int,Int}
    sgp_crpix1::Float64
    sgp_crpix2::Float64
    sgp_lam_scal::Float64
end

"""
    SFD98Map([mapdir])

Schlegel, Finkbeiner and Davis (1998) dust map. 

The first time this is constructed, the data files required will be downloaded 
and stored in a directory following the semantics of 
[DataDeps.jl](https://github.com/oxinabox/datadeps.jl). To avoid being asked to
 download the files, set the environment variable `DATADEPS_ALWAYS_ACCEPT` to 
 `true`. You can also provide the directory of the two requisite files manually 
 instead of relying on DataDeps.jl. Internally, this type keeps the FITS files 
 defining the map open, speeding up repeated queries for E(B-V) values.

# References
[Schlegel, Finkbeiner and Davis (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract)
"""
function SFD98Map(mapdir::String)
    try
        ngp = FITS(joinpath(mapdir, "SFD_dust_4096_ngp.fits"))[1]
        ngp_size = size(ngp)
        ngp_crpix1 = read_key(ngp, "CRPIX1")[1]
        ngp_crpix2 = read_key(ngp, "CRPIX2")[1]
        ngp_lam_scal = read_key(ngp, "LAM_SCAL")[1]
        sgp = FITS(joinpath(mapdir, "SFD_dust_4096_sgp.fits"))[1]
        sgp_size = size(sgp)
        sgp_crpix1 = read_key(sgp, "CRPIX1")[1]
        sgp_crpix2 = read_key(sgp, "CRPIX2")[1]
        sgp_lam_scal = read_key(sgp, "LAM_SCAL")[1]
        SFD98Map(mapdir,
            ngp, ngp_size, ngp_crpix1, ngp_crpix2, ngp_lam_scal,
            sgp, sgp_size, sgp_crpix1, sgp_crpix2, sgp_lam_scal)
    catch
        error("Could not open dust map FITS files in directory $mapdir")
    end
    
end

SFD98Map() = SFD98Map(datadep"sfd98_map")

show(io::IO, map::SFD98Map) = print(io, "SFD98Map($(map.mapdir))")

# Convert from galactic longitude/latitude to lambert pixels.
# See SFD 98 Appendix C. For the 4096x4096 maps, lam_scal = 2048,
# crpix1 = 2048.5, crpix2 = 2048.5.
function galactic_to_lambert(crpix1, crpix2, lam_scal, n, l, b)
    x = lam_scal * sqrt(1 - n * sin(b)) * cos(l) + crpix1
    y = -lam_scal * n * sqrt(1 - n * sin(b)) * sin(l) + crpix2
    return x, y
end

"""
    (dustmap::SFD98Map)(l::Real, b::Real)
    (dustmap::SFD98Map)(l::Quantity, b::Quantity)

Get E(B-V) value from a `SFD98Map` instance at galactic coordinates (`l`, `b`), 
given in radians. Uses bilinear interpolation between pixel values. If `l` and 
`b` are `Unitful.Quantity` they will be converted to radians and the output 
will be given as `UnitfulAstro.mag`. 

# Example

```jldoctest
julia> using DustExtinction

julia> m = SFD98Map();

julia> m(1, 2)
0.013439524544325624

julia> l = 0:0.5:2; b = 0:0.5:2;

julia> m.(l, b)
5-element Array{Float64,1}:
 99.69757461547852    
  0.10180447359074371 
  0.019595484241066132
  0.010238757633890877
  0.01862100327420125 
```

"""
function (dustmap::SFD98Map)(l::Real, b::Real)
    if b >= 0
        hdu = dustmap.ngp
        crpix1 = dustmap.ngp_crpix1
        crpix2 = dustmap.ngp_crpix2
        lam_scal = dustmap.ngp_lam_scal
        xsize, ysize = dustmap.ngp_size
        n = 1
    else
        hdu = dustmap.sgp
        crpix1 = dustmap.sgp_crpix1
        crpix2 = dustmap.sgp_crpix2
        lam_scal = dustmap.sgp_lam_scal
        xsize, ysize = dustmap.sgp_size
        n = -1
    end

    x, y = galactic_to_lambert(crpix1, crpix2, lam_scal, n, l, b)

    # determine integer pixel locations and weights for bilinear interpolation

    xw, xfloor = modf(x)
    x0 = round(Int, xfloor)
    yw, yfloor = modf(y)
    y0 = round(Int, yfloor)

    # handle cases near l = [0, pi/2. pi, 3pi/2] where two pixels
    # are out of bounds. This is made simpler because we know from the
    # galactic_to_lambert() transform that only x or y will be near
    # the image bounds, but not both.
    if x0 == 0
        data = read(hdu, 1, y0:y0 + 1)
        val = (1 - yw) * data[1] + yw * data[2]
    elseif x0 == xsize
        data = read(hdu, xsize, y0:y0 + 1)
        val = (1 - yw) * data[1] + yw * data[2]
    elseif y0 == 0
        data = read(hdu, x0:x0 + 1, 1)
        val = (1 - xw) * data[1] + xw * data[2]
    elseif y0 == ysize
        data = read(hdu, x0:x0 + 1, xsize)
        val = (1 - xw) * data[1] + xw * data[2]
    else
        data = read(hdu, x0:x0 + 1, y0:y0 + 1)
        val = ((1 - xw) * (1 - yw) * data[1, 1] +
               xw       * (1 - yw) * data[2, 1] +
               (1 - xw) * yw       * data[1, 2] +
               xw       * yw       * data[2, 2])
    end
    return val
end

function (dustmap::SFD98Map)(l::Quantity, b::Quantity)
    l_ = ustrip(u"rad", l)
    b_ = ustrip(u"rad", b)
    return dustmap(l_, b_) * u"mag"
end

# Deprecations
@deprecate ebv_galactic(dustmap::SFD98Map, l::Real, b::Real) dustmap(l, b)
@deprecate ebv_galactic(dustmap::SFD98Map, l::AbstractVector, b::AbstractVector) dustmap.(l, b)
