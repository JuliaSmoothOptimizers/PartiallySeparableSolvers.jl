using Documenter
using PartiallySeparableSolvers

using PartiallySeparableSolvers.Mod_partitioned_methods,
  PartiallySeparableSolvers.Mod_TR_CG_part_data
  
makedocs(
  modules = [PartiallySeparableSolvers, Mod_partitioned_methods, Mod_TR_CG_part_data],
  doctest = true,
  # linkcheck = true,
  strict = false,
  format = Documenter.HTML(
    assets = ["assets/style.css"],
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  sitename = "PartiallySeparableSolvers.jl",
  pages = Any["Home" => "index.md", "Tutorial" => "tutorial.md", "Reference" => "reference.md"],
)

deploydocs(repo = "github.com/paraynaud/PartiallySeparableSolvers.jl.git",
  push_preview = true,
  devbranch = "master"
)