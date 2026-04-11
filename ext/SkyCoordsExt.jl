module SkyCoordsExt
    import DustExtinction: SFD98Map, CSFDMap
    using SkyCoords: AbstractSkyCoords, GalCoords, lon, lat

    function (dustmap::Union{SFD98Map, CSFDMap})(s::AbstractSkyCoords)
        s2 = convert(GalCoords, s)
        l, b = lon(s2), lat(s2)
        return dustmap(l, b)
    end

end # module
