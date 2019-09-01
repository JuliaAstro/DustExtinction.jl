using Documenter
using DustExtinction

makedocs(
    modules = [DustExtinction],
    sitename = "DustExtinction.jl",
    doctest = :fix
)

deploydocs(
    repo = "github.com/JuliaAstro/DustExtinction.jl.git",
)
