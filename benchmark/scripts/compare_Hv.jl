using BenchmarkTools
using JuMP, MathOptInterface
using NLPModels, NLPModelsJuMP


using PartiallySeparableSolvers
using PartiallySeparableNLPModel
using CalculusTreeTools
using NLPModels, NLPModelsJuMP


include("generate_problems.jl")

const SUITE = BenchmarkGroup()

problem_collection = create_problems()
global cpt = 1
for (mod, mod_AD) in problem_collection
  println("Hv : itération " * string(cpt) *  "/41")

  n = mod.moi_backend.model_cache.model.num_variables_created

  #= Définition des nlp modèles que nous allons comparer =#
  jump_nlp = NLPModelsJuMP.MathOptNLPModel(mod)
  sps_nlp = PartiallySeparableSolvers.PartionnedNLPModel(jump_nlp)

  #= Définition des points que nous allons tester=#
  x = sps_nlp.meta.x0
  v = ones(eltype(x), length(x))

  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"]["Hv"] = BenchmarkGroup()
  # définition des variables nécessaires

  hv_jump = similar(x)
  hv_sps = similar(x)
  hv_ad = similar(x)
  SUITE["problem $cpt"]["Hv"]["sps"] = @benchmarkable NLPModels.hprod!($sps_nlp, $x, $v, $hv_sps; obj_weight=1.0)
  SUITE["problem $cpt"]["Hv"]["jump"] = @benchmarkable NLPModels.hprod!($jump_nlp, $x, $v, $hv_jump; obj_weight=1.0)
  SUITE["problem $cpt"]["Hv"]["ad"] = @benchmarkable NLPModels.hprod!($mod_AD, $x, $v, $hv_ad; obj_weight=1.0)

  global cpt += 1
end
