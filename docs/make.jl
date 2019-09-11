using Documenter
using DustExtinction

makedocs(
    modules = [DustExtinction],
    sitename = "DustExtinction.jl",
    strict = true
)

deploydocs(
    repo = "github.com/JuliaAstro/DustExtinction.jl.git",
)
