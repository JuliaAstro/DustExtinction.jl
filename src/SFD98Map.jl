using FITSIO

import Base: show

export download_sfd98, SFD98Map, ebv_galactic

# SFD98 Dust Maps

# It would be good to find a more permanent remote location of these maps,
# but they're here for the time being.
const SFD98_BASEURL = "http://sncosmo.github.io/data/dust/"

# Download 4096x4096 E(B-V) maps to $SFD98_DIR directory.
"""
    download_sfd98([destdir])

Download the Schlegel, Finkbeiner and Davis (1998) dust maps to the given
directory. If the directory is ommitted, the `SFD98_DIR` environment variable
is used as the destination directory.
"""
function download_sfd98(destdir::AbstractString)
    for fname in ["SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits"]
        dest = joinpath(destdir, fname)
        if isfile(dest)
            info("$dest already exists, skipping download.")
        else
            download(SFD98_BASEURL * fname, dest)
        end
    end
end

function download_sfd98()
    haskey(ENV, "SFD98_DIR") || error("SFD98_DIR environment variable not set")
    destdir = ENV["SFD98_DIR"]
    download_sfd98(destdir)
end

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

Schlegel, Finkbeiner and Davis (1998) dust map. `mapdir` should be a
directory containing the two FITS files defining the map,
`SFD_dust_4096_[ngp,sgp].fits`. If `mapdir` is omitted, the
`SFD98_DIR` environment variable is used. Internally, this type keeps
the FITS files defining the map open, speeding up repeated queries
for E(B-V) values.
"""
function SFD98Map(mapdir::AbstractString)
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
end

function SFD98Map()
    haskey(ENV, "SFD98_DIR") || error("SFD98_DIR environment variable not set")
    SFD98Map(ENV["SFD98_DIR"])
end

show(io::IO, map::SFD98Map) = print(io, "SFD98Map(\"$(map.mapdir)\")")

# Convert from galactic longitude/latitude to lambert pixels.
# See SFD 98 Appendix C. For the 4096x4096 maps, lam_scal = 2048,
# crpix1 = 2048.5, crpix2 = 2048.5.
function galactic_to_lambert(crpix1, crpix2, lam_scal, n, l, b)
    x = lam_scal * sqrt(1. - n * sin(b)) * cos(l) + crpix1
    y = -lam_scal * n * sqrt(1. - n * sin(b)) * sin(l) + crpix2
    return x, y
end

"""
    ebv_galactic(dustmap::SFD98Map, l::Real, b::Real)
    ebv_galactic(dustmap::SFD98Map, l::Vector{<:Real}, b::Vector{<:Real})

Get E(B-V) value from a `SFD98Map` instance at galactic coordinates
(`l`, `b`), given in radians. `l` and `b` may be Vectors. Uses bilinear
interpolation between pixel values.
"""
function ebv_galactic(dustmap::SFD98Map, l::Real, b::Real)
    if b >= 0.
        hdu = dustmap.ngp
        crpix1 = dustmap.ngp_crpix1
        crpix2 = dustmap.ngp_crpix2
        lam_scal = dustmap.ngp_lam_scal
        xsize, ysize = dustmap.ngp_size
        n = 1.
    else
        hdu = dustmap.sgp
        crpix1 = dustmap.sgp_crpix1
        crpix2 = dustmap.sgp_crpix2
        lam_scal = dustmap.sgp_lam_scal
        xsize, ysize = dustmap.sgp_size
        n = -1.
    end

    x, y = galactic_to_lambert(crpix1, crpix2, lam_scal, n, l, b)

    # determine interger pixel locations and weights for bilinear interpolation
    xfloor = floor(x)
    xw = x - xfloor
    x0 = round(Int, xfloor)
    yfloor = floor(y)
    yw = y - yfloor
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
               xw      * (1 - yw) * data[2, 1] +
               (1 - xw) * yw      * data[1, 2] +
               xw      * yw      * data[2, 2])
    end

    return convert(Float64, val)
end

# Vectorized version
function ebv_galactic(dustmap::SFD98Map, l::Vector{T}, b::Vector{T}) where T <: Real
    m = length(l)
    length(b) == m || error("length of l and b must match")
    result = Array(Float64, m)
    for i = 1:m
        result[i] = ebv_galactic(dustmap, l[i], b[i])
    end
    
    return result
end
