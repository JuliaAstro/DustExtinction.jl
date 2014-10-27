# TODO list:
# 1. Make this optional functionality  - only defined if FITSIO is available
# 2. Function to auto-download maps and place in some standard location?
# 3. Functionality to open files once and access values multiple times.
#    (SFDDustMap type that holds two fits files?)

using FITSIO

# l, b: galactic coordinates in radians
# mapdir: directory where SFD98 dust maps are located
#         with names SFD_dust_4096_[sgp,ngp].fits
function sfd98ebv(l::Real, b::Real, mapdir::String)

    if b >= 0.
        fname = "SFD_dust_4096_ngp.fits"
        sign = 1.
    else
        fname = "SFD_dust_4096_sgp.fits"
        sign = -1.
    end
    fpath = joinpath(mapdir, fname)

    # Open fits file
    f = FITS(fpath, "r")
    
    # Project from galactic longitude/latitude to lambert pixels.
    # (See SFD98).
    npix = mapd.shape[0]        
    x = npix/2 * cos(l) * sqrt(1. - sign*sin(b)) + npix/2 - 0.5
    y = -npix/2 * sign*sin(l) * sqrt(1. - sign*sin(b)) + npix/2 - 0.5

    # cast x, y to int and get the map value there.
    data = f[1][int(x+0.5), int(y+0.5)]
    close(f)
end
