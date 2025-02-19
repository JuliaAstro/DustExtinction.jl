# [Plotting](@id plotting)

DustExtinction.jl provides basic functionality for creating plots with Makie.jl, which can then be composed together to form more complex figures. Below are a few examples using the following primitives from Makie.jl: [`lines`](https://docs.makie.org/stable/reference/plots/lines), [`scatter`](https://docs.makie.org/stable/reference/plots/scatter), and [`heatmap`](https://docs.makie.org/stable/reference/plots/heatmap), along with their mutating equivalents.

!!! note
    By default, all plots adopt the axis limits defined by [`DustExtinction.checkbounds`](@ref), but this is easily modified in the Makie figures.

## Line plots
Makie's line and scatter plots work out-of-the box for all color laws, which can then be used to build up more complex plots.

```@example a
using DustExtinction, CairoMakie
using DustExtinction: aa_to_invum
using Measurements
```

```@example a
wavs = let
    x = range(2_000, 3_000; length=1_000)
    x .± 1e5 * inv.(x)
end # Å
```

```@example a
invum = aa_to_invum.(wavs) # μm
```

```@example a
model = CCM89()
```

```@example a
lines(model)
```

```@example a
lines(model; axis=(; limits=(3, 6, 0, 5)))
```

```@example a
extinction = model.(wavs) # mag
```

```@example a
lines(wavs, extinction)
```

```@example a
wavs_sampled, extinction_sampled = let
    N_samples = 7
    wavs[range(begin, step=end ÷ N_samples; length=N_samples)],
    extinction[range(begin, step=end ÷ N_samples; length=N_samples)]
end
```

```@example a
fig, ax, p = band(wavs, extinction; alpha=0.5, label="model uncertainty")

lines!(ax, wavs, extinction; label="model: CCM89")

scatter!(ax, wavs_sampled, extinction_sampled; color=:orange, label="observations")

# Currently ambiguous for both x and y being Measurements
errorbars!(ax, Measurements.value.(wavs_sampled), extinction_sampled;
    whiskerwidth = 10,
    color = :orange,
    label = "obs. uncertainty",
)

axislegend()
ax.xlabel = "Wavelength [Å]"
ax.ylabel = "A(x) / A(V) [mag]"

fig
```

!!! tip
    See [plotting.jl](https://github.com/JuliaAstro/DustExtinction.jl/blob/docs-makie/docs/plotting.jl) for more plotting examples. The convenience functions defined there are used to generate the other figures shown in this documentation.
