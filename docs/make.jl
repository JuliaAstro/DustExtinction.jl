using Documenter
using DustExtinction

DocMeta.setdocmeta!(DustExtinction, :DocTestSetup, :(using DustExtinction); recursive = true)
include("pages.jl")
makedocs(modules = [DustExtinction],
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
