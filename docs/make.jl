using Documenter
using PartiallySeparableSolvers

makedocs(
  modules = [PartiallySeparableSolvers],
  doctest = true,
  # linkcheck = true,
  strict = true,
  format = Documenter.HTML(
    assets = ["assets/style.css"],
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  sitename = "PartiallySeparableSolvers.jl",
  pages = Any["Home" => "index.md", "Tutorial" => "tutorial.md", "Reference" => "reference.md"],
)

deploydocs(repo = "github.com/paraynaud/PartiallySeparableSolvers.jl.git", devbranch = "main")


# using DocumenterTools
# DocumenterTools.genkeys(user="paraynaud", repo="PartiallySeparableSolvers.jl")