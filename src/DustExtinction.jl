module DustExtinction

using FITSIO  # for SFD98 dust maps

import Base: show

export ccm89,
       od94,
       download_sfd98,
       SFD98Map,
       ebv_galactic


# -----------------------------------------------------------------------------
# Extinction Laws

# Optical coefficients
const ccm89_ca = [1., 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530,
                  0.32999]
const ccm89_cb = [0., 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260,
                  -2.09002]
const od94_ca = [1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647,
                 -0.505]
const od94_cb = [0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805,
                 3.347]

function ccm89like(w::FloatingPoint, r_v, c_a, c_b)
    x = 1.e4 / w
    a = 0.
    b = 0.
    if x < 0.3
        error("wavelength out of range")
    elseif x < 1.1  # Near IR
        y = x^1.61
        a = 0.574y
        b = -0.527y
    elseif x < 3.3  # Optical
        y = x - 1.82
        yn = 1.
        a = c_a[1]
        b = c_b[1]
        for i = 2:length(c_a)
            yn *= y
            a += c_a[i] * yn
            b += c_b[i] * yn
        end
    elseif x < 8.  # UV
        a =  1.752 - 0.316x - (0.104 / ((x-4.67)^2 + 0.341))
        b = -3.090 + 1.825x + (1.206 / ((x-4.62)^2 + 0.263))
        if x > 5.9
            y = x - 5.9
            y2 = y * y
            y3 = y2 * y
            a += -0.04473y2 - 0.009779y3
            b += 0.2130y2 + 0.1207y3
        end
    elseif x < 10.
        y = x - 8.
        y2 = y * y
        y3 = y2 * y
        a = -0.070y3 + 0.137y2 - 0.628y - 1.073
        b = 0.374y3 - 0.420y2 + 4.257y + 13.670
    else
        error("wavelength out of range")
    end

    a + b/r_v
end

ccm89(w::FloatingPoint, r_v) = ccm89like(w, r_v, ccm89_ca, ccm89_cb)
od94(w::FloatingPoint, r_v) = ccm89like(w, r_v, od94_ca, od94_cb)

# Vectorized versions (vectorized on wavelength only)
for f = (:ccm89, :od94)
    @eval begin
        ($f){T<:FloatingPoint}(w::AbstractArray{T,1}, r_v) =
            [ ($f)(w[i], r_v) for i=1:length(w) ]
        ($f){T<:FloatingPoint}(w::AbstractArray{T,2}, r_v) =
            [ ($f)(w[i,j], r_v) for i=1:size(w,1), j=1:size(w,2) ]
        ($f){T<:FloatingPoint}(w::AbstractArray{T}, r_v) =
            reshape([ ($f)(w[i], r_v) for i=1:length(w) ], size(w))
    end
end


# -----------------------------------------------------------------------------
# SFD98 Dust Maps

const SFD98_BASEURL = ("http://www.astro.princeton.edu/~schlegel/dust/" *
                       "dustpub/maps/")

# Download 4096x4096 E(B-V) maps to $SFD98_DIR directory.
function download_sfd98(destdir::String)
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

# SFD98Map: opens both North and South FITS files, reads a few header
# values from each, and keeps them open for future data reads.
type SFD98Map
    mapdir::String
    ngp::ImageHDU
    ngp_size::(Int, Int)
    ngp_crpix1::Float64
    ngp_crpix2::Float64
    ngp_lam_scal::Float64
    sgp::ImageHDU
    sgp_size::(Int, Int)
    sgp_crpix1::Float64
    sgp_crpix2::Float64
    sgp_lam_scal::Float64

    function SFD98Map(mapdir::String)
        ngp = FITS(joinpath(mapdir, "SFD_dust_4096_ngp.fits"))[1]
        ngp_size = size(ngp)
        ngp_crpix1 = readkey(ngp, "CRPIX1")[1]
        ngp_crpix2 = readkey(ngp, "CRPIX2")[1]
        ngp_lam_scal = readkey(ngp, "LAM_SCAL")[1]
        sgp = FITS(joinpath(mapdir, "SFD_dust_4096_sgp.fits"))[1]
        sgp_size = size(sgp)
        sgp_crpix1 = readkey(sgp, "CRPIX1")[1]
        sgp_crpix2 = readkey(sgp, "CRPIX2")[1]
        sgp_lam_scal = readkey(sgp, "LAM_SCAL")[1]
        new(mapdir,
            ngp, ngp_size, ngp_crpix1, ngp_crpix2, ngp_lam_scal,
            sgp, sgp_size, sgp_crpix1, sgp_crpix2, sgp_lam_scal)
    end
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
    x = lam_scal * sqrt(1. - n*sin(b)) * cos(l) + crpix1
    y = -lam_scal * n * sqrt(1. - n*sin(b)) * sin(l) + crpix2
    return x, y
end


# l, b: galactic coordinates in radians
#
# uses bilinear interpolation
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
    x0 = int(xfloor)
    yfloor = floor(y)
    yw = y - yfloor
    y0 = int(yfloor)

    # handle cases near l = [0, pi/2. pi, 3pi/2] where two pixels
    # are out of bounds. This is made simpler because we know from the
    # galactic_to_lambert() transform that only x or y will be near
    # the image bounds, but not both.
    if x0 == 0
        data = hdu[1, y0:y0+1]
        val = (1. - yw) * data[1] + yw * data[2]
    elseif x0 == xsize
        data = hdu[xsize, y0:y0+1]
        val = (1. - yw) * data[1] + yw * data[2]
    elseif y0 == 0
        data = hdu[x0:x0+1, 1]
        val = (1. - xw) * data[1] + xw * data[2]
    elseif y0 == ysize
        data = hdu[x0:x0+1, xsize]
        val = (1. - xw) * data[1] + xw * data[2]
    else
        data = hdu[x0:x0+1, y0:y0+1]
        val = ((1.-xw) * (1.-yw) * data[1, 1] +
               xw      * (1.-yw) * data[2, 1] + 
               (1.-xw) * yw      * data[1, 2] +
               xw      * yw      * data[2, 2])
    end

    return convert(Float64, val)
end

# array version
function ebv_galactic{T <: Real}(dustmap::SFD98Map, l::Vector{T}, b::Vector{T})
    m = length(l)
    length(b) == m || error("length of l and b must match")
    result = Array(Float64, m)
    for i=1:m
        result[i] = ebv_galactic(dustmap, l[i], b[i])
    end
    return result
end

end # module
