# Downloads the SFD98 Map file
const SFD98_BASEURL = "http://sncosmo.github.io/data/dust/"

destdir = haskey(ENV, "SFD98_DIR") ? ENV["SFD98_DIR"] : @__DIR__

try
    for fname in ["SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits"]
        dest = joinpath(destdir, fname)
        if isfile(dest)
            info("$dest already exists, skipping download.")
        else
            download(SFD98_BASEURL * fname, dest)
        end
    end
catch
    @warn "There was an error downloading the SFD98 Dust Maps. Please re-build or manually download the maps and set the \"SFD98_DIR\" environment variable."
end
