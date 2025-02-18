module MakieExt
    using DustExtinction
    import DustExtinction: bounds, aa_to_invum, ExtinctionLaw
    import Makie

    function Makie.convert_arguments(P::Makie.PointBased, law::ExtinctionLaw)
        aa = range(bounds(law)...; length=1_000)
        m = map(law, aa)
        invum = map(aa_to_invum, aa)
        return Makie.convert_arguments(P, invum, m) # (Plot type, x, y)
    end

    Makie.plottype(::ExtinctionLaw) = Makie.Lines
end
