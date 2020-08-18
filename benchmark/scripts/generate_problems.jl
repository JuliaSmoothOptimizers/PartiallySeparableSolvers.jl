using OptimizationProblems, JuMP, ADNLPModels



#=
liste = arwhead, bdqrtic, brybnd, chainwoo, chnrosnb_mod, cosine , cragglvy, curly, curly10, curly20, curly30, dixmaane, dixmaanf, dixmaang, dixmaan,
dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp, dixon3dq, dqdrtic, edensch, eg2, engval1, errinros_mod, extrosnb,
freuroth, genhumps, liarwhd, morebv, noncvxu2, noncvxun, nondia, nondquar, palmer1d, palmer1c, palmer2c, palmer3c, power, quartc, sbrybnd, schmvett, scosine,
sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, tridia, vardim


liste_a_tester = clplatea, clplateb, clplatec, fminsrf2

liste_sps_elem = arglina, arglinb, arglinc

liste_sps_manque_op = NZF1 (log), broydn7d (abs)

pas Partiellement séparable = power
=#

include("define_ADNLPModel.jl")

function create_problems(n :: Int)
  problem_collection = Vector{Tuple{JuMP.Model,ADNLPModels.RADNLPModel}}(undef,0)
  push!(problem_collection, (OptimizationProblems.arwhead(n), arwhead_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.bdqrtic(n), bdqrtic_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.brybnd(n), brybnd_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.chainwoo(n - (n %4)), chainwoo_ADNLPModel(n - (n %4))) )
  push!(problem_collection, (OptimizationProblems.chnrosnb_mod(n), chnrosnb_mod_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.cragglvy(n), cragglvy_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.curly(n), curly_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.curly10(n), curly10_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.curly20(n), curly20_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.curly30(n), curly30_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.dixmaane(n - (n %3)), dixmaane_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaanf(n - (n %3)), dixmaanf_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaang(n - (n %3)), dixmaang_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaanh(n - (n %3)), dixmaanh_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaani(n - (n %3)), dixmaani_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaanj(n - (n %3)), dixmaanj_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaank(n - (n %3)), dixmaank_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaanl(n - (n %3)), dixmaanl_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaanm(n - (n %3)), dixmaanm_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaann(n - (n %3)), dixmaann_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaano(n - (n %3)), dixmaano_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixmaanp(n - (n %3)), dixmaanp_ADNLPModel(n - (n %3))) )
  push!(problem_collection, (OptimizationProblems.dixon3dq(n), dixon3dq_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.dqdrtic(n), dqdrtic_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.dqrtic(n), dqrtic_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.edensch(n), edensch_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.eg2(n), eg2_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.errinros_mod(n), errinros_mod_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.extrosnb(n), extrosnb_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.freuroth(n), freuroth_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.genhumps(n), genhumps_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.liarwhd(n), liarwhd_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.morebv(n), morebv_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.noncvxu2(n), noncvxu2_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.noncvxun(n), noncvxun_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.nondia(n), nondia_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.nondquar(n), nondquar_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.quartc(n), quartc_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.sbrybnd(n), sbrybnd_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.tridia(n), tridia_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.vardim(n), vardim_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.scosine(n), scosine_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.sinquad(n), sinquad_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.vardim(n), vardim_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.sparsine(n), sparsine_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.sparsqur(n), sparsqur_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.srosenbr(n), srosenbr_ADNLPModel(n)) )
  push!(problem_collection, (OptimizationProblems.woods(n - (n %4)), woods_ADNLPModel(n - (n %4))) )








  # push!(problem_collection, OptimizationProblems.arwhead(n))
  # push!(problem_collection, OptimizationProblems.bdqrtic(n))
  # push!(problem_collection, OptimizationProblems.brybnd(n))
  # push!(problem_collection, OptimizationProblems.chainwoo(n))
  # push!(problem_collection, OptimizationProblems.chnrosnb_mod(n))
  # push!(problem_collection, OptimizationProblems.cosine(n))
  # push!(problem_collection, OptimizationProblems.cragglvy(n))
  # push!(problem_collection, OptimizationProblems.curly(n))
  # push!(problem_collection, OptimizationProblems.curly10(n))
  # push!(problem_collection, OptimizationProblems.curly20(n))
  # push!(problem_collection, OptimizationProblems.curly30(n))
  # push!(problem_collection, OptimizationProblems.dixmaane(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaanf(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaang(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaanh(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaani(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaanj(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaank(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaanl(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaanm(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaann(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaano(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixmaanp(n - (n %3)))
  # push!(problem_collection, OptimizationProblems.dixon3dq(n))
  # push!(problem_collection, OptimizationProblems.dqdrtic(n))
  # push!(problem_collection, OptimizationProblems.qdrtic(n))
  # push!(problem_collection, OptimizationProblems.edensch(n))
  # push!(problem_collection, OptimizationProblems.eg2(n))
  # push!(problem_collection, OptimizationProblems.engval1(n))
  # push!(problem_collection, OptimizationProblems.errinros_mod(n))
  # push!(problem_collection, OptimizationProblems.extrosnb(n))
  # push!(problem_collection, OptimizationProblems.freuroth(n))
  # push!(problem_collection, OptimizationProblems.genhumps(n))
  # push!(problem_collection, OptimizationProblems.liarwhd(n))
  # push!(problem_collection, OptimizationProblems.morebv(n))
  # push!(problem_collection, OptimizationProblems.noncvxu2(n))
  # push!(problem_collection, OptimizationProblems.noncvxun(n))
  # push!(problem_collection, OptimizationProblems.nondia(n))
  # push!(problem_collection, OptimizationProblems.nondquar(n))
  # push!(problem_collection, OptimizationProblems.quartc(n))
  # push!(problem_collection, OptimizationProblems.sbrybnd(n))
  # push!(problem_collection, OptimizationProblems.tridia(n))
  # push!(problem_collection, OptimizationProblems.vardim(n))
  # push!(problem_collection, OptimizationProblems.scosine(n))
  # push!(problem_collection, OptimizationProblems.sinquad(n))
  # push!(problem_collection, OptimizationProblems.sparsine(n))
  # push!(problem_collection, OptimizationProblems.sparsqur(n))
  # push!(problem_collection, OptimizationProblems.srosenbr(n))
  # push!(problem_collection, OptimizationProblems.woods(n))


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

n = 10
@inline create_problems() = create_problems(n)
