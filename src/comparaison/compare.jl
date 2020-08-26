repo_model = "models/"
include(repo_model * "chained_wood.jl")
include( repo_model * "rosenbrock.jl")
include( repo_model * "chained_powel.jl")
include( repo_model * "chained_cragg_levy.jl")
include( repo_model * "generalisation_Brown.jl")



using JSOSolvers, SolverBenchmark, SolverTools, Plots, Printf, DataFrames, NLPModels
using ADNLPModels
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

"""
    second_criteria(n, time, obj, neval_grad, neval_Hv)
Returns a kind of score =time / ( n_obj + 5*neval_grad + 5*neval_Hv). This criteria evluate the efficiency of the base function used in a solver
"""
second_criteria(n, time, n_obj, n_grad, n_hess) = time / (n_obj + 5*n_grad + 5*n_hess)
function my_criteria_time(d :: DataFrames.DataFrame)
    nvar = d.nvar
    iter = d.iter
    n_eval_obj = d.neval_obj
    n_eval_grad = d.neval_grad
    n_eval_hprod = d.neval_hprod
    elapsed_time = d.elapsed_time
    my_criteria = []
    for i in 1:length(nvar)
      push!(my_criteria, second_criteria(nvar[i], elapsed_time[i], n_eval_obj[i], n_eval_grad[i], n_eval_hprod[i]))
    end
    d.time_sur_obj_5grad_5Hv = my_criteria
end


"""
    second_criteria(n_obj, neval_grad)
Returns the percentage of step keep by our solver
"""
third_criteria(n_obj, n_grad) = 1/((n_grad/n_obj)*100)
function acceptance_criteria(d :: DataFrames.DataFrame)
  n_eval_obj = d.neval_obj
  n_eval_grad = d.neval_grad
  my_criteria = []
  for i in 1:length(n_eval_obj)
    push!(my_criteria, third_criteria(n_eval_obj[i], n_eval_grad[i]))
  end
  d.inverse_pourcentage_pas_accepte = my_criteria
end




println(" \n\n génération des problemes")
# n_array = [100,500,1000,2000,5000,10000,25000,50000,100000]
# n_array = [1000,2000,5000, 10000]
# problems = create_all_problems(n_array)
# problems2 = create_all_problems2(n_array)
include("../../benchmark/scripts/generate_problems.jl")
problems = create_JuMP_models()
problems2 = create_ADNLP_models()


println("\n\ndéfinition des solver\n\n")



const atol = 1.0e-5
const rtol = 1.0e-6
const max_time = 600.0
const max_eval = 1500


solver2 = Dict{Symbol,Function}(
  :trunk_adnlpmodel => (prob; kwargs...) -> JSOSolvers.trunk(prob; kwargs...),
  :lsr1_adnlpmodel => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LSR1Model(prob);kwargs...),
  :lbfgs_adnlpmodel => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LBFGSModel(prob);kwargs...)
)


solver = Dict{Symbol,Function}(
:trunk => ((prob;kwargs...) -> JSOSolvers.trunk(prob;kwargs...)),
:trunk_lsr1 => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LSR1Model(prob); kwargs...),
:trunk_lbfgs => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LBFGSModel(prob); kwargs...),
:trunk_SPS => (prob; kwargs...) -> JSOSolvers.trunk(PartiallySeparableSolvers.PartionnedNLPModel(prob); kwargs...),
# :my_lbfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LBFGS(prob;kwargs...)),
# :my_lsr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LSR1(prob;kwargs...)),
:p_bfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.PBFGS(prob; kwargs...)),
:p_sr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.PSR1(prob; kwargs...)),
:p_bs => ((prob;kwargs...) -> PartiallySeparableSolvers.PBS(prob; kwargs...)),
# :p_trunk => ((prob;kwargs...) -> PartiallySeparableSolvers.PTRUNK(prob; kwargs...))
)




keys_hess = [:trunk, :trunk_adnlpmodel, :trunk_SPS]
keys_bfgs = [:trunk_lbfgs, :lbfgs_adnlpmodel, :p_bfgs]
keys_sr1 =  [:trunk_lsr1, :lsr1_adnlpmodel, :p_sr1, :p_bs ]



#= Lancement du benchmark sur les problèmes générés, sur les solvers défini dans la variable solvers =#
println("lancement des benchmmarks ADNLPModel")
#lancement de bmark_solver sur les ADNLPModels
stats2 = bmark_solvers(solver2, problems2; max_time=max_time, max_eval = max_eval, atol=atol, rtol=rtol)
#récupération des clés
keys_stats2 = keys(stats2)

println("lancement des benchmmarks NLPModelJuMP")
#lancement de bmark_solver sur les NLPModelJUMP
stats = bmark_solvers(solver, problems; max_time=max_time, max_eval = max_eval, atol=atol, rtol=rtol)

#on ajoute les Dataframes de stats2 à stats
for i in keys_stats2
  stats[i] = stats2[i]
end






println("affichage du profile des solvers par rapport au problèmes")
# création des colonnes liés aux critères que j'ai défini
keys_stats = keys(stats)
for i in keys_stats
  my_criteria_timeless(stats[i])
  my_criteria_time(stats[i])
  acceptance_criteria(stats[i])
end

#Construction de tables particulières pour chaque classe de solvers.
stats_hess = Dict{Symbol,DataFrame}([])
for i in keys_hess
  stats_hess[i] = stats[i]
end

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
selected_fields = [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :obj_5grad_5Hv, :time_sur_obj_5grad_5Hv, :inverse_pourcentage_pas_accepte]
for i in keys_stats
  println(stdout, "\n\n\n" * string(i) )
  markdown_table(stdout, stats[i], cols=selected_fields)
end

#= Ecriture des résultats dans un fichier au format markdown=#
println("écriture des résultats markdown")
location_md = string("src/comparaison/results/result_bench_md.txt")
io = open(location_md,"w")
close(io)
io = open(location_md,"w+")

for i in keys_stats
  println(io, "\n\n\n" * string(i) )
  markdown_table(io, stats[i], cols=selected_fields)
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
  latex_table(io, stats[i], cols=selected_fields)
end
close(io)



println("ecriture des profiles")
#pas d'affichage des profils, server amdahl
ENV["GKSwstype"]=100

println("écriture de tous les profiles")
p_iter = SolverBenchmark.performance_profile(stats, df -> df.iter)
savefig(p_iter, "src/comparaison/results/profiles/iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats, df -> df.elapsed_time)
savefig(p_time, "src/comparaison/results/profiles/time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, "src/comparaison/results/profiles/obj_5grad_5Hv.pdf")
p_snd_crit = SolverBenchmark.performance_profile(stats, df -> df.time_sur_obj_5grad_5Hv)
savefig(p_snd_crit, "src/comparaison/results/profiles/time_sur_obj_5grad_5Hv.pdf")
p_thd_crit = SolverBenchmark.performance_profile(stats, df -> df.inverse_pourcentage_pas_accepte)
savefig(p_thd_crit, "src/comparaison/results/profiles/inverse_pourcentage_pas_accepte.pdf")

println("écriture des profiles hess like")

repo_hess = "src/comparaison/results/profiles/hess_like/"
p_iter = SolverBenchmark.performance_profile(stats_hess, df -> df.iter)
savefig(p_iter, repo_hess * "iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats_hess, df -> df.elapsed_time)
savefig(p_time, repo_hess * "time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats_hess, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, repo_hess * "obj_5grad_5Hv.pdf")
p_snd_crit = SolverBenchmark.performance_profile(stats_hess, df -> df.time_sur_obj_5grad_5Hv)
savefig(p_snd_crit, repo_hess * "time_sur_obj_5grad_5Hv.pdf")
p_thd_crit = SolverBenchmark.performance_profile(stats_hess, df -> df.inverse_pourcentage_pas_accepte)
savefig(p_thd_crit, repo_hess * "inverse_pourcentage_pas_accepte.pdf")


println("écriture des profiles BFGS like")

repo_bfgs = "src/comparaison/results/profiles/bfgs_like/"
p_iter = SolverBenchmark.performance_profile(stats_bfgs, df -> df.iter)
savefig(p_iter, repo_bfgs * "iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats_bfgs, df -> df.elapsed_time)
savefig(p_time, repo_bfgs * "time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats_bfgs, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, repo_bfgs * "obj_5grad_5Hv.pdf")
p_snd_crit = SolverBenchmark.performance_profile(stats_bfgs, df -> df.time_sur_obj_5grad_5Hv)
savefig(p_snd_crit, repo_bfgs * "time_sur_obj_5grad_5Hv.pdf")
p_thd_crit = SolverBenchmark.performance_profile(stats_bfgs, df -> df.inverse_pourcentage_pas_accepte)
savefig(p_thd_crit, repo_bfgs * "inverse_pourcentage_pas_accepte.pdf")


println("écriture des profiles SR1 like")

repo_sr1 = "src/comparaison/results/profiles/sr1_like/"
p_iter = SolverBenchmark.performance_profile(stats_sr1, df -> df.iter)
savefig(p_iter, repo_sr1 * "iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats_sr1, df -> df.elapsed_time)
savefig(p_time, repo_sr1 * "time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats_sr1, df -> df.obj_5grad_5Hv)
savefig(p_fst_crit, repo_sr1 * "obj_5grad_5Hv.pdf")
p_snd_crit = SolverBenchmark.performance_profile(stats_sr1, df -> df.time_sur_obj_5grad_5Hv)
savefig(p_snd_crit, repo_sr1 * "time_sur_obj_5grad_5Hv.pdf")
p_thd_crit = SolverBenchmark.performance_profile(stats_sr1, df -> df.inverse_pourcentage_pas_accepte)
savefig(p_thd_crit, repo_sr1 * "inverse_pourcentage_pas_accepte.pdf")


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
# files of the repo models. n Is the size of the problems. Trhis function generates RADNLPModels.
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
