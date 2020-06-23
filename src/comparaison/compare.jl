repo_model = "models/"
include(repo_model * "chained_wood.jl")
include( repo_model * "rosenbrock.jl")
include( repo_model * "chained_powel.jl")
include( repo_model * "chained_cragg_levy.jl")
include( repo_model * "generalisation_Brown.jl")



using JSOSolvers, SolverBenchmark, SolverTools, Plots, Printf, DataFrames

using PartiallySeparableSolvers



"""
    create_all_problems(n)
Function that create the problem that I will use with bmark_solver. I use function to generate automaticaly some Models, theses are defined in other
files in the same repo. n Is the size of the problems.
"""
function create_all_problems(nb_var_array :: Vector{Int})
  problem_array = []
  for i in nb_var_array
    (m_ros, evaluator,obj) = create_Rosenbrock_JuMP_Model(i)
    (m_chained, evaluator,obj) = create_chained_wood_JuMP_Model(i)
    (m_powel, evaluator,obj) = create_chained_Powel_JuMP_Model(i)
    (m_cragg_levy, evaluator,obj) = create_chained_cragg_levy_JuMP_Model(i)
    (m_brown, evaluator,obj) = create_generalisation_Brown(i)
    push!(problem_array, MathOptNLPModel(m_ros, name="Ros "*string(i)), MathOptNLPModel(m_chained, name="ChainedWood "*string(i)), MathOptNLPModel(m_powel, name="ChainedPowel "*string(i)), MathOptNLPModel(m_cragg_levy, name="CraggLevy "*string(i)), MathOptNLPModel(m_brown, name="Brown "*string(i)))
    # init = create_initial_point_Rosenbrock(i)
    # push!(problem_array, (MathOptNLPModel(m_ros), init))
  end
  return problem_array
end


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
    d.first_criteria = my_criteria
end

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
    d.second_criteria = my_criteria
end


third_criteria(n_obj, n_grad) = (n_grad/n_obj)*100
function acceptance_criteria(d :: DataFrames.DataFrame)
  n_eval_obj = d.neval_obj
  n_eval_grad = d.neval_grad
  my_criteria = []
  for i in 1:length(n_eval_obj)
    push!(my_criteria, third_criteria(n_eval_obj[i], n_eval_grad[i]))
  end
  d.third_criteria = my_criteria
end


# stats[:trunk].champ = []
# markdown_table(stdout, stats[:trunk], cols=[champ])




println(" \n\n génération des problemes")
n_array = [100,500,1000,2000,5000,10000,20000,40000,60000,80000,100000]
# n_array = [10,20,30]
# n_array = [1000,2000]
# n_array = [100,200]
# n_array = [100,500,1000,2000,5000]
problems = create_all_problems(n_array)


println("\n\ndéfinition des solver\n\n")


solver = Dict{Symbol,Function}(
  :trunk => ((prob;kwargs...) -> JSOSolvers.trunk(prob;kwargs...)),
  :trunk_lsr1 => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LSR1Model(prob); kwargs...),
  :my_lbfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LBFGS(prob;kwargs...)),
  :my_lsr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LSR1(prob;kwargs...)),
  :p_bfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.PBFGS(prob; kwargs...)),
  :p_sr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.PSR1(prob; kwargs...)),
  :p_bs => ((prob;kwargs...) -> PartiallySeparableSolvers.PBS(prob; kwargs...))
)

const atol = 1.0e-5
const rtol = 1.0e-6
const max_time = 300.0


#= Lancement du benchmark sur les problèmes générés, sur les solvers défini dans la variable solvers =#
println("lancement des benchmmarks")

stats = bmark_solvers(solver, problems; max_time=max_time, max_eval = 5000, atol=atol, rtol=rtol)

println("affichage du profile des solvers par rapport au problèmes")
# performance_profile(stats, df->df.elapsed_time)
# performance_profile(stats, df->df.iter)
keys_stats = keys(stats)
for i in keys_stats
  my_criteria_timeless(stats[i])
  my_criteria_time(stats[i])
  acceptance_criteria(stats[i])
end

println("affichage des tables")
markdown_table(stdout, stats[:trunk], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
markdown_table(stdout, stats[:trunk_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
markdown_table(stdout, stats[:my_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
markdown_table(stdout, stats[:my_lbfgs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
markdown_table(stdout, stats[:p_bfgs], cols=[ :name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
markdown_table(stdout, stats[:p_sr1], cols=[  :name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
markdown_table(stdout, stats[:p_bs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])


# stats[:trunk].champ = []
# markdown_table(stdout, stats[:trunk], cols=[champ])

# error("fin")
#= Ecriture des résultats dans un fichier au format markdown=#
println("écriture des résultats markdown")
location_md = string("src/comparaison/results/result_bench_md.txt")
io = open(location_md,"w")
close(io)
io = open(location_md,"w+")
println(io, "\n\n\nTrunk" )
markdown_table(io, stats[:trunk], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nTrunk LSR1" )
markdown_table(io, stats[:trunk_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nMon LSR1" )
markdown_table(io, stats[:my_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nMon LBFGS" )
markdown_table(io, stats[:my_lbfgs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nPBFGS" )
markdown_table(io, stats[:p_bfgs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nPSR1" )
markdown_table(io, stats[:p_sr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nPBS" )
markdown_table(io, stats[:p_bs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])

close(io)


#= Ecriture des résultats dans un fichier au format latex=#
println("écriture des résultats latex")
location_latex = string("src/comparaison/results/result_bench_latex.txt")
io = open(location_latex,"w")
close(io)
io = open(location_latex,"w+")
println(io, "\n\n\nTrunk" )
latex_table(io, stats[:trunk], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nTrunk LSR1" )
latex_table(io, stats[:trunk_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nMon LSR1" )
latex_table(io, stats[:my_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nMon LBFGS" )
latex_table(io, stats[:my_lbfgs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nPBFGS" )
latex_table(io, stats[:p_bfgs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nPSR1" )
latex_table(io, stats[:p_sr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
println(io, "\n\n\nPBS" )
latex_table(io, stats[:p_bs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad, :neval_hprod, :first_criteria, :second_criteria, :third_criteria])
close(io)
