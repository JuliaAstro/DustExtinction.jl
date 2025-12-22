module MakieExt
    using DustExtinction
    import DustExtinction: bounds, aa_to_invum, ExtinctionLaw
    import Makie as M

    # Ex: plot(CCM89())
    function M.convert_arguments(P::M.PointBased, law::ExtinctionLaw)
        aa = range(bounds(law)...; length=1_000)
        m = map(law, aa)
        invum = map(aa_to_invum, aa)
        return M.convert_arguments(P, invum, m)
    end
    M.plottype(::ExtinctionLaw) = M.Lines

    # Ex: plot(wavs, CCM89())
    function M.convert_arguments(P::M.PointBased, x, law::ExtinctionLaw)
        m = map(law, x)
        M.convert_arguments(P, x, m)
    end
    M.plottype(x, ::ExtinctionLaw) = M.Lines

    # Ex: heatmap(SFD98Map())
    function M.convert_arguments(P::M.CellGrid, dustmap::SFD98Map)
        l = range(-3, 3; length=400)
        b = range(-1, 1; length=300)
        m = [dustmap(li, bj) for li in l, bj in b]
        return M.convert_arguments(P, l, b, m)
    end
    M.plottype(::SFD98Map) = M.Heatmap

    # Ex: heatmap(lrange, brange, SFD98Map())
    function M.convert_arguments(P::M.CellGrid, lrange, brange, dustmap::SFD98Map)
        m = [dustmap(li, bj) for li in lrange, bj in brange]
        return M.convert_arguments(P, lrange, brange, m)
    end
    M.plottype(x, y, ::SFD98Map) = M.Heatmap
end
