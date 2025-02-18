# [Plotting](@id plotting)

DustExtinction.jl provides basic functionality for creating plots with Makie.jl, which can then be composed together to form more complex figures. Below are a few examples using the following primitives from Makie.jl: [`lines`](https://docs.makie.org/stable/reference/plots/lines), [`scatter`](https://docs.makie.org/stable/reference/plots/scatter), and [`heatmap`](https://docs.makie.org/stable/reference/plots/heatmap), along with their mutating equivalents.

!!! note
    By default, all plots adopt the axis limits defined by [`DustExtinction.checkbounds`](@ref), but this is easily modified in the Makie figures.

## Line plots
Makie's line and scatter plots work out-of-the box for all color laws, which can then be used to build up more complex plots.

```@example plot
using DustExtinction, CairoMakie
fig, ax1, p1 = lines(CCM89())
ax2, p2 = scatter(fig[1, 2], F04(); color=:orange, markersize=5)
linkaxes!(ax1, ax2)
xlims!(3, 8)
fig
```

!!! tip
    See [plotting.jl](https://github.com/JuliaAstro/DustExtinction.jl/blob/docs-makie/docs/plotting.jl) for more plotting examples. These functions are used to generate the other figures shown in this documentation.
