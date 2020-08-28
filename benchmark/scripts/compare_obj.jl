using BenchmarkTools
using NLPModels, NLPModelsJuMP


using PartiallySeparableSolvers

include("generate_problems.jl")

const SUITE = BenchmarkGroup()
SUITE["jump"] = BenchmarkGroup()
SUITE["ad"] = BenchmarkGroup()
SUITE["sps"] = BenchmarkGroup()

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


  SUITE["jump"]["problem $cpt"] = BenchmarkGroup()
  SUITE["ad"]["problem $cpt"] = BenchmarkGroup()
  SUITE["sps"]["problem $cpt"] = BenchmarkGroup()

  #calcul de la fonction objectif
  println("\t sps")
  SUITE["sps"]["problem $cpt"] = @benchmarkable NLPModels.obj($sps_nlp, $x)
  println("\t JuMP")
  SUITE["jump"]["problem $cpt"] = @benchmarkable NLPModels.obj($jump_nlp, $x)
  println("\t AD")
  SUITE["ad"]["problem $cpt"] = @benchmarkable NLPModels.obj($mod_AD, $x)


  global cpt += 1
end
