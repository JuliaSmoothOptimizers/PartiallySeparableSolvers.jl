using JSOSolvers, SolverBenchmark, SolverTools, Plots, Printf, DataFrames, NLPModels
using ADNLPModels
using ProfileView
using PartiallySeparableSolvers



"""
    first_criteria(n, obj, neval_grad, neval_Hv)
Returns a kind of score = n_obj + 5*neval_grad + 5*neval_Hv. This criteria is a approximation of the cost of the function
"""
first_criteria(n, n_obj, n_grad, n_hess) = n_obj + 5*n_grad + 5*n_hess
function my_criteria_timeless(d :: DataFrames.DataFrame)
    nvar = d.nvar
    n_eval_obj = d.neval_obj
    n_eval_grad = d.neval_grad
    n_eval_hprod = d.neval_hprod
    my_criteria = []
    for i in 1:length(nvar)
      push!(my_criteria, first_criteria(nvar[i], n_eval_obj[i], n_eval_grad[i], n_eval_hprod[i]))
    end
    d.obj_5grad_5Hv = my_criteria
end




println(" \n\n génération des problemes")
current_path = pwd()
include(current_path * "/benchmark/scripts/generate_problems.jl")


n_array = [100,200,400,800,1600,3200,6400,12800]
mapreduce((n -> create_JuMP_models(n)), vcat, n_array) # generate the problem for the sizes of n_array
# problems = create_JuMP_models()
# problems2 = create_ADNLP_models()


println("\n\ndéfinition des solver\n\n")



const atol = 1.0e-5
const rtol = 1.0e-6
const max_time = 300.0
const max_eval = 100000


# solver2 = Dict{Symbol,Function}(
  # :trunk_Hv_adnlpmodel => (prob; kwargs...) -> JSOSolvers.trunk(prob; kwargs...),
  # :trunk_lsr1_adnlpmodel => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LSR1Model(prob);kwargs...),
  # :trunk_lbfgs_adnlpmodel => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LBFGSModel(prob);kwargs...)
# )


solver = Dict{Symbol,Function}(
# :trunk_Hv_JuMP => ((prob;kwargs...) -> JSOSolvers.trunk(prob;kwargs...)),
# :trunk_lsr1_JuMP => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LSR1Model(prob); kwargs...),
# :trunk_lbfgs_JuMP => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LBFGSModel(prob); kwargs...),
# :trunk_Hv_SPS => (prob; kwargs...) -> JSOSolvers.trunk(PartiallySeparableSolvers.PartionnedNLPModel(prob); kwargs...),
:my_lbfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LBFGS(prob;kwargs...)),
:my_lsr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LSR1(prob;kwargs...)),
:bfgs_SPS => ((prob;kwargs...) -> PartiallySeparableSolvers.PBFGS(prob; kwargs...)),
:sr1_SPS => ((prob;kwargs...) -> PartiallySeparableSolvers.PSR1(prob; kwargs...)),
:bs_SPS => ((prob;kwargs...) -> PartiallySeparableSolvers.PBS(prob; kwargs...)),
# :p_trunk => ((prob;kwargs...) -> PartiallySeparableSolvers.PTRUNK(prob; kwargs...))
)




keys_hess = [:trunk_Hv_JuMP, :trunk_Hv_adnlpmodel, :trunk_Hv_SPS]
# keys_bfgs = [:trunk_lbfgs_JuMP, :bfgs_SPS, :bs_SPS]
keys_bfgs = [:my_lbfgs, :bfgs_SPS, :bs_SPS]
# keys_sr1 =  [:trunk_lsr1_JuMP, :sr1_SPS, :bs_SPS ]
keys_sr1 =  [:my_lsr1, :sr1_SPS, :bs_SPS ]


#= Lancement du benchmark sur les problèmes générés, sur les solvers défini dans la variable solvers =#

println("lancement des benchmmarks NLPModelJuMP")
#lancement de bmark_solver sur les NLPModelJUMP
stats = bmark_solvers(solver, problems; max_time=max_time, max_eval = max_eval, atol=atol, rtol=rtol)

# println("lancement des benchmmarks ADNLPModel")
# #lancement de bmark_solver sur les ADNLPModels
# stats2 = bmark_solvers(solver2, problems2; max_time=max_time, max_eval = max_eval, atol=atol, rtol=rtol)

# #récupération des clés
# keys_stats2 = keys(stats2)

# #on ajoute les Dataframes de stats2 à stats
# for i in keys_stats2
#   stats[i] = stats2[i]
# end






println("affichage du profile des solvers par rapport au problèmes")
# création des colonnes liés aux critères que j'ai défini
keys_stats = keys(stats)
for i in keys_stats
  my_criteria_timeless(stats[i])
end

#Construction de tables particulières pour chaque classe de solvers.

stats_bfgs = Dict{Symbol,DataFrame}([])
for i in keys_bfgs
  stats_bfgs[i] = stats[i]
end

stats_sr1 = Dict{Symbol,DataFrame}([])
for i in keys_sr1
  stats_sr1[i] = stats[i]
end

println("affichage des tables")
#selection des champs à affichier
selected_fields = [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :obj_5grad_5Hv]
for i in keys_stats
  println(stdout, "\n\n\n" * string(i) )
  pretty_stats(stdout, stats[i][!, [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod]], tf=markdown)
end

#= Ecriture des résultats dans un fichier au format markdown=#
println("écriture des résultats markdown")
location_md = string("src/comparaison/results/result_bench_md.txt")
io = open(location_md,"w")
close(io)
io = open(location_md,"w+")

for i in keys_stats
  println(io, "\n\n\n" * string(i) )
  pretty_stats(io, stats[i][!, [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod]], tf=markdown)
end
close(io)



#= Ecriture des résultats dans un fichier au format latex=#
println("écriture des résultats latex")
location_latex = string("src/comparaison/results/result_bench_latex.txt")
io = open(location_latex,"w")
close(io)
io = open(location_latex,"w+")
for i in keys_stats
  println(io, "\n\n\n" * string(i) )
  pretty_latex_stats(io, stats[i][!, [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod]])
end
close(io)



println("ecriture des profiles")
#pas d'affichage des profils, server amdahl
ENV["GKSwstype"]=100

println("écriture de tous les profiles")
repo_global = "src/comparaison/results/amdahl/profile/"
p_iter = SolverBenchmark.performance_profile(stats, df -> df.iter)
savefig(p_iter, repo_global*"iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats, df -> df.elapsed_time)
savefig(p_time, repo_global*"time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, repo_global*"obj_5grad_5Hv.pdf")



println("écriture des profiles BFGS like")

repo_bfgs = "src/comparaison/results/amdahl/bfgs_like/"
p_iter = SolverBenchmark.performance_profile(stats_bfgs, df -> df.iter)
savefig(p_iter, repo_bfgs * "iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats_bfgs, df -> df.elapsed_time)
savefig(p_time, repo_bfgs * "time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats_bfgs, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, repo_bfgs * "obj_5grad_5Hv.pdf")

println("écriture des profiles SR1 like")

repo_sr1 = "src/comparaison/results/amdahl/sr1_like/"
p_iter = SolverBenchmark.performance_profile(stats_sr1, df -> df.iter)
savefig(p_iter, repo_sr1 * "iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats_sr1, df -> df.elapsed_time)
savefig(p_time, repo_sr1 * "time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats_sr1, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, repo_sr1 * "obj_5grad_5Hv.pdf")

println("Fin des écritures")




#=
ancien Code
=#

#
# """
#     create_all_problems(n)
# Function that create the problem that I will use with bmark_solver. I use function to generate automaticaly some Models, theses are defined in other
# files of the repo models. n Is the size of the problems. Trhis function generates NLPMoodelsJuMP.
# """
# function create_all_problems(nb_var_array :: Vector{Int})
#   problem_array = []
#   for i in nb_var_array
#     push!(problem_array, create_rosenbrock_JuMPModel(i))
#     push!(problem_array, create_chained_wood_JuMPModel(i))
#     push!(problem_array, create_chained_powel_JuMPModel(i))
#     push!(problem_array, create_cragg_levy_JuMPModel(i))
#     push!(problem_array, create_generalisation_brown_JuMPModel(i))
#   end
#   return problem_array
# end
#
#
# """
#     create_all_problems2(n)
# Function that create the problem that I will use with bmark_solver. I use function to generate automaticaly some Models, theses are defined in other
# files of the repo models. n Is the size of the problems. Trhis function generates ADNLPModels.
# """
# function create_all_problems2(nb_var_array :: Vector{Int})
#   problem_array = []
#   for i in nb_var_array
#     push!(problem_array, create_rosenbrock_ADNLPModel(i))
#     push!(problem_array, create_chained_wood_ADNLPModel(i))
#     push!(problem_array, create_chained_powel_ADNLPModel(i))
#     push!(problem_array, create_cragg_levy_ADNLPModel(i))
#     push!(problem_array, create_generalisation_brown_ADNLPModel(i))
#   end
#   return problem_array
# end
# keys_hess = [:trunk, :trunk_adnlpmodel, :trunk_SPS, :p_trunk]
# keys_bfgs = [ :lbfgs_adnlpmodel, :my_lbfgs, :p_bfgs]
# keys_sr1 =  [:trunk_lsr1, :lsr1_adnlpmodel, :my_lsr1, :p_sr1, :p_bs ]
