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
  println("obj : itération " * string(cpt) *  "/41")

  n = mod.moi_backend.model_cache.model.num_variables_created


  #= Définition des nlp modèles que nous allons comparer =#
  println("\t Model JuMP")
  jump_nlp = NLPModelsJuMP.MathOptNLPModel(mod)
  println("\t Model Partitonned")
  sps_nlp = PartiallySeparableSolvers.PartionnedNLPModel(jump_nlp)

  #= Définition des points que nous allons tester=#
  x = sps_nlp.meta.x0
  v = ones(eltype(x), length(x))

  SUITE["problem $cpt"] = BenchmarkGroup()
  SUITE["problem $cpt"] = BenchmarkGroup()
  # SUITE["problem $cpt"]["obj"] = BenchmarkGroup()
  SUITE["problem $cpt"]["jump"] = BenchmarkGroup()
  SUITE["problem $cpt"]["ad"] = BenchmarkGroup()
  SUITE["problem $cpt"]["sps"] = BenchmarkGroup()

  # définition des variables nécessaires

  #calcul de la fonction objectif
  println("\t sps")
  SUITE["problem $cpt"]["sps"] = @benchmarkable NLPModels.obj($sps_nlp, $x)
  println("\t JuMP")
  SUITE["problem $cpt"]["jump"] = @benchmarkable NLPModels.obj($jump_nlp, $x)
  println("\t AD")
  SUITE["problem $cpt"]["ad"] = @benchmarkable NLPModels.obj($mod_AD, $x)

  # SUITE["problem $cpt"]["obj"]["sps"] = @benchmarkable NLPModels.obj($sps_nlp, $x)
  # SUITE["problem $cpt"]["obj"]["ad"] = @benchmarkable NLPModels.obj($sps_nlp, $x)
  # SUITE["problem $cpt"]["obj"]["jump"] = @benchmarkable NLPModels.obj($sps_nlp, $x)

  global cpt += 1
end
