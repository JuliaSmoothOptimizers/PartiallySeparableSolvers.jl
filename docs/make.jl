using Documenter
using PartiallySeparableSolvers

using PartiallySeparableSolvers.ModPartitionedMethods,
  PartiallySeparableSolvers.ModTrustRegionPartitionedData,
  PartiallySeparableSolvers.Mod_PQN,
  PartiallySeparableSolvers.Mod_ab_partitioned_data

makedocs(
  modules = [
    PartiallySeparableSolvers,
    ModPartitionedMethods,
    ModTrustRegionPartitionedData,
    Mod_PQN,
    Mod_ab_partitioned_data,
  ],
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
