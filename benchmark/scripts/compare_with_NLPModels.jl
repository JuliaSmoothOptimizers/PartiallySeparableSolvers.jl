using BenchmarkTools
using JuMP, MathOptInterface
using NLPModels, NLPModelsJuMP

using PartiallySeparableSolvers

include("generate_problems.jl")

const SUITE = BenchmarkGroup()

problem_collection = create_problems()
global cpt = 1
for (mod, mod_AD) in problem_collection
  n = mod.moi_backend.model_cache.model.num_variables_created


  #= Définition des nlp modèles que nous allons comparer =#
  jump_nlp = NLPModelsJuMP.MathOptNLPModel(mod)
  sps_nlp = PartiallySeparableSolvers.PartionnedNLPModel(jump_nlp)

  #= Définition des points que nous allons tester=#
  x = sps_nlp.meta.x0
  v = ones(eltype(x), length(x))

  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"]["obj"] = BenchmarkGroup()
  SUITE["problem $cpt"]["grad"] = BenchmarkGroup()
  SUITE["problem $cpt"]["Hv"] = BenchmarkGroup()
  # définition des variables nécessaires

  #calcul de la fonction objectif
  SUITE["problem $cpt"]["obj"]["sps"] = @benchmarkable NLPModels.obj($sps_nlp, $x)
  SUITE["problem $cpt"]["obj"]["jump"] = @benchmarkable NLPModels.obj($jump_nlp, $x)
  SUITE["problem $cpt"]["obj"]["ad"] = @benchmarkable NLPModels.obj($mod_AD, $x)

  #calcul du gradient sous format gradient élémentaire
  grad_sps = similar(x)
  grad_jump = similar(x)
  grad_AD = similar(x)
  SUITE["problem $cpt"]["grad"]["sps"] = @benchmarkable NLPModels.grad!($sps_nlp, $x, $grad_sps)
  SUITE["problem $cpt"]["grad"]["jump"] = @benchmarkable NLPModels.grad!($jump_nlp, $x, $grad_jump)
  SUITE["problem $cpt"]["grad"]["ad"] = @benchmarkable NLPModels.grad!($mod_AD, $x, $grad_AD)


  hv_jump = similar(x)
  hv_sps = similar(x)
  hv_AD = similar(x)
  SUITE["problem $cpt"]["Hv"]["sps"] = @benchmarkable NLPModels.hprod!($sps_nlp, $x, $v, $hv_sps; obj_weight=1.0)
  SUITE["problem $cpt"]["Hv"]["jump"] = @benchmarkable NLPModels.hprod!($jump_nlp, $x, $v, $hv_jump; obj_weight=1.0)
  SUITE["problem $cpt"]["Hv"]["ad"] = @benchmarkable NLPModels.hprod!($mod_AD, $x, $v, $hv_AD; obj_weight=1.0)

  global cpt += 1
end
