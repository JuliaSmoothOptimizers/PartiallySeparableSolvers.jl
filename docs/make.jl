using Documenter
using PartiallySeparableSolvers

using PartiallySeparableSolvers.ModPartitionedMethods,
  PartiallySeparableSolvers.ModTrustRegionPartitionedData

makedocs(
  modules = [PartiallySeparableSolvers, ModPartitionedMethods, ModTrustRegionPartitionedData],
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

deploydocs(
  repo = "github.com/paraynaud/PartiallySeparableSolvers.jl.git",
  push_preview = true,
  devbranch = "master",
)
