# [Plotting](@id plotting)

DustExtinction.jl supports Makie.jl plots and provides the following recipes for visualizing color laws and dust maps.

### Building blocks
Makie's line and scatter plots work out-of-the box for all color laws, which can then be used to build up more complex plots.

```@example
using DustExtinction, CairoMakie

fig, ax1, p1 = lines(CCM89())
ax2, p2 = scatter(fig[1, 2], F04(); color=:orange, markersize=5)
linkaxes!(ax1, ax2)
fig
```

### Parameter average models
For more complete plots that include automatic axis labels and legends, use [`lplot`](@ref).

```@example
using DustExtinction, CairoMakie

lplot(CCM89())
```

```@example
using DustExtinction, CairoMakie

lplot(CCM89)
```

### Shape models

### Mixture models

### Dust maps

```@example
using DustExtinction
using CairoMakie

dplot()
```

See [`Color laws`] or [`Dust Maps`] for more usage examples.

## API/Reference

```@docs
lplot
dplot
```
