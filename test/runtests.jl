using DustExtinction
using Test
using Documenter

include("ccm89.jl")
include("cal00.jl")
include("sfd98.jl")
doctest(DustExtinction)

println("Tests passed.")
