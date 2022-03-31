using JSOSolvers, SolverBenchmark, SolverTools, Plots, Printf, DataFrames, NLPModels
using ADNLPModels
using PartiallySeparableSolvers


"""
    first_criteria(n, obj, neval_grad, neval_Hv)
Returns a kind of score = n_obj + 5*neval_grad + 5*neval_Hv. This criteria is a approximation of the cost of the function
"""
first_criteria(n, n_obj, n_grad, n_hess) = Float64(n_obj + 5*n_grad + 5*n_hess)
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
		@show my_criteria
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
# n_array = [1000,2000,5000, 10000]
# problems = create_all_problems(n_array)
current_path = pwd()
include(current_path * "/benchmark/scripts/generate_problems.jl")
n=1000
problems = create_JuMP_models(n) #MathOptNLPModel

println("\n\ndéfinition des solver\n\n")
const atol = 1.0e-5
const rtol = 1.0e-6
const max_time = 600.0
const max_eval = 5000


solver = Dict{Symbol,Function}(
:PBFGS => ((prob;kwargs...) -> PartiallySeparableSolvers.PUS(prob; name=:pbfgs, kwargs...)),
:PLBFGS => ((prob;kwargs...) -> PartiallySeparableSolvers.PUS(prob; name=:plbfgs, kwargs...)),
:PSR1 => ((prob;kwargs...) -> PartiallySeparableSolvers.PUS(prob; name=:psr1, kwargs...)),
:PLSR1 => ((prob;kwargs...) -> PartiallySeparableSolvers.PUS(prob; name=:plsr1, kwargs...)),
:PSE => ((prob;kwargs...) -> PartiallySeparableSolvers.PUS(prob; name=:pse, kwargs...)),
:PLSE => ((prob;kwargs...) -> PartiallySeparableSolvers.PUS(prob; name=:plse, kwargs...)),
)

solver2 = Dict{Symbol,Function}(
:LBFGS => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LBFGS(prob;kwargs...)),
)


# keys_hess = [:trunk_Hv_JuMP, :trunk_Hv_adnlpmodel, :trunk_Hv_SPS]
keys_solvers = [:LBFGS, :PBFGS, :PLBFGS, :PSR1, :PLSR1, :PSE, :PLSE]


#= Lancement du benchmark sur les problèmes générés, sur les solvers défini dans la variable solvers =#

println("lancement des benchmmarks NLPModelJuMP")
#lancement de bmark_solver sur les NLPModelJUMP
stats = bmark_solvers(solver, problems; max_time=max_time, max_eval = max_eval, atol=atol, rtol=rtol, max_iter=max_eval)

println("lancement des benchmmarks ADNLPModel")
#lancement de bmark_solver sur les ADNLPModels
stats2 = bmark_solvers(solver2, problems; max_time=max_time, max_eval = max_eval, atol=atol, rtol=rtol)

#récupération des clés
keys_stats2 = keys(stats2)

#on ajoute les Dataframes de stats2 à stats
for i in keys_stats2
  stats[i] = stats2[i]
end

println("affichage du profile des solvers par rapport au problèmes")
# création des colonnes liés aux critères que j'ai défini
keys_stats = keys(stats)
for i in keys_stats
	println(i)
  my_criteria_timeless(stats[i])
  my_criteria_time(stats[i])
  acceptance_criteria(stats[i])
end

#Construction de tables particulières pour chaque classe de solvers.

# stats_keys = Dict{Symbol,DataFrame}([])
# for i in keys_solvers
#   stats_keys[i] = stats[i]
# end
# cost(df) = (df.status .!= :success) * Inf + df.t

println("affichage des tables")
#selection des champs à affichier
selected_fields = [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :obj_5grad_5Hv, :time_sur_obj_5grad_5Hv, :inverse_pourcentage_pas_accepte]
for i in keys_stats
  println(stdout, "\n\n\n" * string(i) )
  pretty_stats(stdout, stats[i][!, [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod]], tf=tf_markdown)
end

#= Ecriture des résultats dans un fichier au format markdown=#
println("écriture des résultats markdown")
location_md = string("src/comparaison/results/plse/result_bench.md")
io = open(location_md,"w")
close(io)
io = open(location_md,"w+")
for i in keys_stats
  println(io, "\n\n\n" * string(i) )
  pretty_stats(io, stats[i][!, [:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod]], tf=tf_markdown)
end
close(io)


#= Ecriture des résultats dans un fichier au format latex=#
println("écriture des résultats latex")
location_latex = string("src/comparaison/results/plse/result_bench_latex.tex")
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
p_iter = SolverBenchmark.performance_profile(stats, df -> (df.status .== :max_eval) * Inf + df.iter; legend=:bottomright )
savefig(p_iter, "src/comparaison/results/plse/profiles/iter_profile.pdf")
p_time = SolverBenchmark.performance_profile(stats, df -> (df.status .== :max_eval) * Inf + df.elapsed_time; legend=:bottomright )
savefig(p_time, "src/comparaison/results/plse/profiles/time_profile.pdf")
p_fst_crit = SolverBenchmark.performance_profile(stats, df -> (df.status .== :max_eval) * Inf + df.obj_5grad_5Hv; legend=:bottomright )
savefig(p_fst_crit, "src/comparaison/results/plse/profiles/obj_5grad_5Hv.pdf")
# p_snd_crit = SolverBenchmark.performance_profile(stats, df -> df.time_sur_obj_5grad_5Hv; legend=:bottomright )
# savefig(p_snd_crit, "src/comparaison/results/plse/profiles/time_sur_obj_5grad_5Hv.pdf")
# p_thd_crit = SolverBenchmark.performance_profile(stats, df -> df.inverse_pourcentage_pas_accepte; legend=:bottomright )
# savefig(p_thd_crit, "src/comparaison/results/plse/profiles/inverse_pourcentage_pas_accepte.pdf")


println("Fin des écritures")