module SkyCoordsExt
    import DustExtinction: SFD98Map
    using SkyCoords: AbstractSkyCoords, GalCoords, lon, lat

    function (dustmap::SFD98Map)(s::AbstractSkyCoords)
        s2 = convert(GalCoords, s)
        l, b = lon(s2), lat(s2)
        return dustmap(l, b)
    end

end # module