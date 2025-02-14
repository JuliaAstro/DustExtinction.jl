module MakieExt
    using DustExtinction
    import DustExtinction: bounds, aa_to_invum, ExtinctionLaw, lplot, lplot!, mplot, mplot!, dplot, dplot!
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

    Rv law plot with automatic axis labels.
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

    Rv law series plot with automatic axis labels.
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
        lplot(law::FM90; args...)

    Fittable law plot with automatic axis labels.
    """
    lplot(law::FM90; args...) = lines(law;
        axis = (;
            xlabel = rich("x [μm", superscript("-1"), "]"),
            ylabel = rich("(E(λ) - V) / E(B - V)"),
        ),
        args...,
    )

    """
        lplot(law::FM90; args...)

    Fittable law series plot with automatic axis labels.
    """
    function lplot(law::Type{<:FM90}; args...)
        # m1
        fig, ax, p = lplot(law(); label="total")

        # m2
        lines!(ax, law(c3=0.0, c4=0.0); label="linear term")

        # m3
        lines!(ax, law(c1=0.0, c2=0.0, c4=0.0); label="bump term")

        # m4
        lines!(ax, law(c1=0.0, c2=0.0, c3=0.0); label="FUV rise term")

        # for Rv in (2.0, 3.1, 4.0, 5.0, 6.0)
        #   lines!(ax, law(Rv); label=rich("Rv = $(Rv)"), args...)
        # end

        axislegend(ax; position=:lt)

        fig
    end

    """
        mplot(law::Type{G16}, Rvs, f_A::Real; args...)

    Mixture law plot with automatic axis labels.
    """
    function mplot(law::Type{G16}, Rvs, f_A::Real; args...)
        # Dummy plot
        fig, ax, p = lplot(law())

        for Rv in Rvs
            lines!(ax, law(; Rv, f_A); label=rich("Rv = $(Rv)"), args...)
        end

        axislegend(ax, "f_A = $(f_A)"; position=:lt)

        fig
    end

    """
        mplot(law::Type{G16}, Rvs, f_A::Real; args...)

    Mixture law series plot with automatic axis labels.
    """
    function mplot(law::Type{G16}, Rv::Real, f_As; args...)
        # Dummy plot
        fig, ax, p = lplot(law())

        for f_A in f_As
            lines!(ax, law(; Rv, f_A); label=rich("f_A = $(f_A)"), args...)
        end

        axislegend(ax, "Rv = $(Rv)"; position=:lt)

        fig
    end

    """
        dplot(dustmap=SFD98Map(); lrange=(-3, 3), brange=(-1, 1))

    Heatmap of the given `dustmap`, with galactic longitude ``(l)`` on the x-axis and galactic latitude ``(b)`` on the y-axis. Angles are displayed in degrees by default. The plot axis ranges for both are set by `lrange` and `brange`, respectively.
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
