using Documenter
using DustExtinction

makedocs(
    modules = [DustExtinction],
    sitename = "DustExtinction.jl",
)

deploydocs(
    repo = "github.com/JuliaAstro/DustExtinction.jl.git",
)
