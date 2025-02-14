module MakieExt
    using DustExtinction
    import DustExtinction: bounds, aa_to_invum, ExtinctionLaw, lplot, lplot!, dplot, dplot!
    using Makie: Makie, convert_arguments, PointBased, lines, lines!, heatmap, heatmap!, Colorbar, rich, superscript, axislegend

    function Makie.convert_arguments(P::PointBased, law::ExtinctionLaw)
        aa = range(bounds(law)...; length=1_000)
        m = map(law, aa)
        invum = map(aa_to_invum, aa)
        return convert_arguments(P, invum, m) # (Plot type, x, y)
    end

    Makie.plottype(::ExtinctionLaw) = Lines

    """
        lplot(law::ExtinctionLaw; args...)

    Color law plot with automatic axis labels.
    """
    lplot(law::ExtinctionLaw; args...) = lines(
        law;
        axis = (;
            xlabel = rich("x [μm", superscript("-1"), "]"),
            ylabel = rich("A(x) / A(V)"),
        ),
        args...,
    )

    """
        lplot(law::Type{<:ExtinctionLaw}; args...)

    Color law series plot with automatic axis labels.
    """
    function lplot(law::Type{<:ExtinctionLaw}; args...)
        # Dummy plot
        fig, ax, p = lplot(law())

        for Rv in (2.0, 3.1, 4.0, 5.0, 6.0)
            lines!(ax, law(Rv); label=rich("Rv = $(Rv)"), args...)
        end

        axislegend(ax; position=:lt)

        fig
    end

    """
        dplot(dustmap=SFD98Map(); lrange=(-3, 3), brange=(-1, 1))

    Plot a heatmap of the given `dustmap`, with galactic longitude ``(l)`` on the x-axis and galactic latitude ``(b)`` on the y-axis. Angles are displayed in degrees by default. The plot axis ranges for both are set by `lrange` and `brange`, respectively.
    """
    function dplot(dustmap=SFD98Map(); lrange=(-3, 3), brange=(-1, 1))
        l = range(lrange..., length=400)
        b = range(brange..., length=300)
        m = [dustmap(li, bj) for li in l, bj in b]

        fig, ax, p = heatmap(l, b, m; colorrange=(0, 3), colormap=:cividis)
        ax.xlabel = "l (°)"
        ax.ylabel = "b (°)"
        Colorbar(fig[1, 2], p; label="E(B - V)")

        fig
    end
end
