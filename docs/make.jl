using Documenter
using DustExtinction
using CairoMakie
using SkyCoords

CairoMakie.activate!(type="png", px_per_unit=3)

DocMeta.setdocmeta!(DustExtinction, :DocTestSetup, :(using DustExtinction); recursive = true)
MakieExt = Base.get_extension(DustExtinction, :MakieExt)
SkyCoordsExt = Base.get_extension(DustExtinction, :SkyCoordsExt)

makedocs(;
    modules = [DustExtinction, MakieExt, SkyCoordsExt],
    sitename = "DustExtinction.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = [
            "assets/favicon.ico",
        ],
        canonical = "https://JuliaAstro.org/DustExtinction/stable/",
    ),
    authors = "Kyle Barbary, Mosé Giordano, Miles Lucas",
    pages = [
        "Home" => "index.md",
        "Color Laws" => "color_laws.md",
        "Dust Maps" => "dust_maps.md",
        "Plotting" => "plotting.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/JuliaAstro/DustExtinction.jl.git",
    push_preview = true,
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
