module MakieExt
    using DustExtinction
    import DustExtinction: bounds, aa_to_invum, ExtinctionLaw
    import Makie

    # Ex: plot(CCM89())
    function Makie.convert_arguments(P::Makie.PointBased, law::ExtinctionLaw)
        aa = range(bounds(law)...; length=1_000)
        m = map(law, aa)
        invum = map(aa_to_invum, aa)
        return Makie.convert_arguments(P, invum, m)
    end
    Makie.plottype(::ExtinctionLaw) = Makie.Lines

    # Ex: plot(wavs, CCM89())
    function Makie.convert_arguments(P::Makie.PointBased, x, law::ExtinctionLaw)
        m = map(law, x)
        Makie.convert_arguments(P, x, m)
    end
    Makie.plottype(x, ::ExtinctionLaw) = Makie.Lines
end
