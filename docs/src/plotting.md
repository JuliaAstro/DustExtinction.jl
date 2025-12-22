# [Plotting](@id plotting)

DustExtinction.jl is designed to work automatically with many of Makie.jl's core plotting functions, which can then be composed together to form more complex figures. Below we show a few practical applications.

!!! note
    By default, all plots adopt the wavelength limits (in units of 1/μm) defined by [`DustExtinction.bounds`](@ref).

## Model plot example
For a given [`DustExtinction.ExtinctionLaw`](@ref), Makie's usual PointBased plotting functions (e.g., [lines](https://docs.makie.org/stable/reference/plots/lines), [scatter](https://docs.makie.org/stable/reference/plots/scatter), [stairs](https://docs.makie.org/stable/reference/plots/stairs), etc.) should work right out-of-the-box:

```@example a
using DustExtinction, CairoMakie
```

```@example a
model = CCM89()
```

```@example a
# Automatic limits defined by model
lines(model; linewidth=3, color=:cornflowerblue) # Or plot(model; [plot options])
```

```@example a
# Custom limits
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2]; limits=(3, 6, 0, 5))

lines!(ax1, model)
stairs!(ax2, model; color=:orange)

linkyaxes!(ax1, ax2)

fig
```

A predefined vector of wavelengths can also be passed to these plotting functions directly. Let's combine this with Makie's integration with Measurements.jl to visualize the underlying uncertainty in our data:

```@example a
using Measurements, Unitful, UnitfulAstro

wavs = let
    x = range(2_000, 3_000; length=1_000)
    x = x .± 1e5 * inv.(ustrip.(x)) # "Å"
end

extinction = model.(wavs) # mag

wavs_sampled, extinction_sampled = let
    N_samples = 7
    wavs[range(begin, step=end ÷ N_samples; length=N_samples)],
    extinction[range(begin, step=end ÷ N_samples; length=N_samples)]
end

fig, ax, p = band(wavs, extinction; alpha=0.5, label="model uncertainty")

lines!(ax, wavs, extinction; label="model: CCM89")

scatter!(ax, wavs_sampled, extinction_sampled; color=:orange, label="observations")

# Currently ambiguous for both x and y being Measurements
# so we focus on the y-uncertainty instead
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

## Dust map example
A [`heatmap`](https://docs.makie.org/stable/reference/plots/heatmap#heatmap) plot also works automatically for [`DustExtinction.SFD98Map`](@ref):

```@example a
dustmap = SFD98Map()

heatmap(dustmap; colorrange=(0, 3), colormap=:cividis) # Or plot(dustmap; [plot kwargs])
```

Similarly to the extinction law plots, we can create our own custom dust map plots. Here is an example using [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/)'s integration:

```@example a
using Unitful

lrange = range(-1, 1; length=400)u"°"
brange = range(-0.25, 0.25; length=300)u"°"

fig, ax, p = heatmap(lrange, brange, dustmap; colorrange=(0, 3), colormap=:cividis)

ax.xlabel = "l [°]"
ax.ylabel = "b [°]"

Colorbar(fig[1, 2], p; label="E(B - V) [mag]")

fig
```

!!! tip
    See [plotting.jl](https://github.com/JuliaAstro/DustExtinction.jl/blob/docs-makie/docs/src/plotting.jl) for more plotting examples. The convenience functions defined there are used to generate the other figures shown in this documentation.
