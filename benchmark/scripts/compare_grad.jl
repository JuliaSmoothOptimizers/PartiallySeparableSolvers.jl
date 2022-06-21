using BenchmarkTools
using JuMP, MathOptInterface
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
  println("grad : itération " * string(cpt) * "/41")

  n = mod.moi_backend.model_cache.model.num_variables_created

  #= Définition des nlp modèles que nous allons comparer =#
  jump_nlp = NLPModelsJuMP.MathOptNLPModel(mod)
  sps_nlp = PartiallySeparableSolvers.PartionnedNLPModel(jump_nlp)

  #= Définition des points que nous allons tester=#
  x = sps_nlp.meta.x0
  v = ones(eltype(x), length(x))

  grad_sps = similar(x)
  grad_jump = similar(x)
  grad_ad = similar(x)

  SUITE["jump"]["problem $cpt"] = BenchmarkGroup()
  SUITE["ad"]["problem $cpt"] = BenchmarkGroup()
  SUITE["sps"]["problem $cpt"] = BenchmarkGroup()

  #calcul de la fonction objectif
  println("\t sps")
  SUITE["sps"]["problem $cpt"] = @benchmarkable NLPModels.grad!($sps_nlp, $x, $grad_sps)
  println("\t JuMP")
  SUITE["jump"]["problem $cpt"] = @benchmarkable NLPModels.grad!($jump_nlp, $x, $grad_jump)
  println("\t AD")
  SUITE["ad"]["problem $cpt"] = @benchmarkable NLPModels.grad!($mod_AD, $x, $grad_ad)

  global cpt += 1
end
