using Documenter
using DustExtinction

makedocs(
    modules = [DustExtinction],
    sitename = "DustExtinction.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
        assets = [
            "assets/front_image.png",
        ]),
    authors = "Kyle Barbary, Mosé Giordano, Miles Lucas",
    strict = true,
    pages = [
        "Home" => "index.md",
        "Color Laws" => "color_laws.md",
        "Dust Maps" => "dust_maps.md"
    ],
)

############################################
# Set up for pushing preview docs from PRs #
############################################
if haskey(ENV, "TRAVIS_PULL_REQUEST") && ENV["TRAVIS_PULL_REQUEST"] != "false"
    @info "Pushing preview docs."
    PR = ENV["TRAVIS_PULL_REQUEST"]
    # Overwrite Documenter's function for generating the versions.js file
    foreach(Base.delete_method, methods(Documenter.Writers.HTMLWriter.generate_version_file))
    Documenter.Writers.HTMLWriter.generate_version_file(_, _) = nothing
    # Overwrite necessary environment variables to trick Documenter to deploy
    ENV["TRAVIS_PULL_REQUEST"] = "false"
    ENV["TRAVIS_BRANCH"] = "master"
    deploydocs(
        devurl = "preview-PR$(PR)",
        repo = "github.com/JuliaAstro/DustExtinction.jl.git",
    )
    exit(0)
end

deploydocs(repo = "github.com/JuliaAstro/DustExtinction.jl.git")
