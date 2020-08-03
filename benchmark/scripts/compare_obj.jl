using BenchmarkTools
using JuMP, MathOptInterface
using NLPModels, NLPModelsJuMP


using PartiallySeparableSolvers
using PartiallySeparableNLPModel
using CalculusTreeTools
using NLPModels, NLPModelsJuMP

using OptimizationProblems

#=
liste = arwhead, bdqrtic, brybnd, chainwoo, chnrosnb_mod, cosine , cragglvy, curly, curly10, curly20, curly30, dixmaane, dixmaanf, dixmaang, dixmaan,
dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp, dixon3dq, dqdrtic, edensch, eg2, engval1, errinros_mod, extrosnb,
freuroth, genhumps, liarwhd, morebv, noncvxu2, noncvxun, nondia, nondquar, palmer1d, palmer1c, palmer2c, palmer3c, power, quartc, sbrybnd, schmvett, scosine,
sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, tridia, vardim


liste_a_tester = clplatea, clplateb, clplatec, fminsrf2

liste_sps_elem = arglina, arglinb, arglinc

liste_sps_manque_op = NZF1 (log), broydn7d (abs)
=#



function create_problems(n :: Int)
  problem_collection = Vector{JuMP.Model}(undef,0)
  push!(problem_collection, OptimizationProblems.arwhead(n))
  push!(problem_collection, OptimizationProblems.bdqrtic(n))
  push!(problem_collection, OptimizationProblems.brybnd(n))
  push!(problem_collection, OptimizationProblems.chainwoo(n))
  push!(problem_collection, OptimizationProblems.chnrosnb_mod(n))
  push!(problem_collection, OptimizationProblems.cosine(n))
  push!(problem_collection, OptimizationProblems.cragglvy(n))
  push!(problem_collection, OptimizationProblems.curly(n))
  push!(problem_collection, OptimizationProblems.curly10(n))
  push!(problem_collection, OptimizationProblems.curly20(n))
  push!(problem_collection, OptimizationProblems.curly30(n))
  push!(problem_collection, OptimizationProblems.dixmaane(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaanf(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaang(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaanh(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaani(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaanj(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaank(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaanl(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaanm(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaann(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixmaanp(n - (n %3)))
  push!(problem_collection, OptimizationProblems.dixon3dq(n))
  push!(problem_collection, OptimizationProblems.dqdrtic(n))
  push!(problem_collection, OptimizationProblems.edensch(n))
  push!(problem_collection, OptimizationProblems.eg2(n))
  push!(problem_collection, OptimizationProblems.engval1(n))
  push!(problem_collection, OptimizationProblems.errinros_mod(n))
  push!(problem_collection, OptimizationProblems.extrosnb(n))
  push!(problem_collection, OptimizationProblems.freuroth(n))
  push!(problem_collection, OptimizationProblems.genhumps(n))
  push!(problem_collection, OptimizationProblems.liarwhd(n))
  push!(problem_collection, OptimizationProblems.morebv(n))
  push!(problem_collection, OptimizationProblems.noncvxu2(n))
  push!(problem_collection, OptimizationProblems.noncvxun(n))
  push!(problem_collection, OptimizationProblems.nondia(n))
  push!(problem_collection, OptimizationProblems.nondquar(n))
  push!(problem_collection, OptimizationProblems.power(n))
  push!(problem_collection, OptimizationProblems.quartc(n))
  push!(problem_collection, OptimizationProblems.sbrybnd(n))
  push!(problem_collection, OptimizationProblems.tridia(n))
  push!(problem_collection, OptimizationProblems.vardim(n))
  push!(problem_collection, OptimizationProblems.scosine(n))
  push!(problem_collection, OptimizationProblems.sinquad(n))
  push!(problem_collection, OptimizationProblems.sparsine(n))
  push!(problem_collection, OptimizationProblems.sparsqur(n))
  push!(problem_collection, OptimizationProblems.srosenbr(n))

# problème using opératuer fac, réglé mais non push du rpertoire
  # push!(problem_collection, OptimizationProblems.tointgss(n))
  # push!(problem_collection, OptimizationProblems.tquartic(n))

  #problèmes algorithmique
  # push!(problem_collection, OptimizationProblems.palmer1d(n))
  # push!(problem_collection, OptimizationProblems.palmer1c(n))
  # push!(problem_collection, OptimizationProblems.palmer2c(n))
  # push!(problem_collection, OptimizationProblems.palmer3c(n))
  # push!(problem_collection, OptimizationProblems.schmvett(n))


  return problem_collection
end

const SUITE = BenchmarkGroup()

problem_collection = create_problems(100)
global cpt = 1
for mod in problem_collection
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

  # définition des variables nécessaires

  #calcul de la fonction objectif
  SUITE["problem $cpt"]["obj"]["sps"] = @benchmarkable NLPModels.obj($sps_nlp, $x)
  SUITE["problem $cpt"]["obj"]["jump"] = @benchmarkable NLPModels.obj($jump_nlp, $x)


  global cpt += 1
end
