using Documenter
using DustExtinction
using CairoMakie

CairoMakie.activate!(type="png", px_per_unit=3)

DocMeta.setdocmeta!(DustExtinction, :DocTestSetup, :(using DustExtinction); recursive = true)
MakieExt = Base.get_extension(DustExtinction, :MakieExt)
include("pages.jl")
makedocs(modules = [DustExtinction, MakieExt],
    sitename = "DustExtinction.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "master",
    ),
    authors = "Kyle Barbary, Mosé Giordano, Miles Lucas",
    pages = pages,
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/JuliaAstro/DustExtinction.jl.git",
    push_preview = true,
    branch = "gh-pages",
    devbranch = "master",
)
