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
  println("Hv : itération " * string(cpt) *  "/41")

  n = mod.moi_backend.model_cache.model.num_variables_created

  #= Définition des nlp modèles que nous allons comparer =#
  jump_nlp = NLPModelsJuMP.MathOptNLPModel(mod)
  sps_nlp = PartiallySeparableSolvers.PartionnedNLPModel(jump_nlp)

  #= Définition des points que nous allons tester=#
  x = sps_nlp.meta.x0
  v = ones(eltype(x), length(x))

  hv_sps = similar(x)
  hv_jump = similar(x)
  hv_ad = similar(x)

  SUITE["jump"]["problem $cpt"] = BenchmarkGroup()
  SUITE["ad"]["problem $cpt"] = BenchmarkGroup()
  SUITE["sps"]["problem $cpt"] = BenchmarkGroup()

  #calcul de la fonction objectif
  println("\t sps")
  SUITE["sps"]["problem $cpt"] = @benchmarkable NLPModels.hprod!($sps_nlp, $x, $v, $hv_sps; obj_weight=1.0)
  println("\t JuMP")
  SUITE["jump"]["problem $cpt"] = @benchmarkable NLPModels.hprod!($jump_nlp, $x, $v, $hv_jump; obj_weight=1.0)
  println("\t AD")
  SUITE["ad"]["problem $cpt"] = @benchmarkable NLPModels.hprod!($mod_AD, $x, $v, $hv_ad; obj_weight=1.0)

  global cpt += 1
end
