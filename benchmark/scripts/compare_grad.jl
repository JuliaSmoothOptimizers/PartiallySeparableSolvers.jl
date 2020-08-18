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
  n = mod.moi_backend.model_cache.model.num_variables_created


  #= Définition des nlp modèles que nous allons comparer =#
  jump_nlp = NLPModelsJuMP.MathOptNLPModel(mod)
  sps_nlp = PartiallySeparableSolvers.PartionnedNLPModel(jump_nlp)

  #= Définition des points que nous allons tester=#
  x = sps_nlp.meta.x0
  v = ones(eltype(x), length(x))

  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"]["grad"] = BenchmarkGroup()
  # définition des variables nécessaires

    #calcul du gradient sous format gradient élémentaire
  grad_sps = similar(x)
  grad_jump = similar(x)
  grad_ad = similar(x)
  SUITE["problem $cpt"]["grad"]["sps"] = @benchmarkable NLPModels.grad!($sps_nlp, $x, $grad_sps)
  SUITE["problem $cpt"]["grad"]["jump"] = @benchmarkable NLPModels.grad!($jump_nlp, $x, $grad_jump)
  SUITE["problem $cpt"]["grad"]["ad"] = @benchmarkable NLPModels.grad!($mod_AD, $x, $grad_ad)

  global cpt += 1
end
