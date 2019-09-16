using Documenter
using DustExtinction

makedocs(modules = [DustExtinction],
    sitename = "DustExtinction.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Kyle Barbary, Mosé Giordano, Miles Lucas",
    strict = true,
    pages = [
        "Home" => "index.md",
        "Color Laws" => "color_laws.md",
        "Dust Maps" => "dust_maps.md"
    ],
)

deploydocs(repo = "github.com/JuliaAstro/DustExtinction.jl.git")
