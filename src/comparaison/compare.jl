repo_model = "models/"
include(repo_model * "chained_wood.jl")
include( repo_model * "rosenbrock.jl")
include( repo_model * "chained_powel.jl")
include( repo_model * "chained_cragg_levy.jl")
include( repo_model * "generalisation_Brown.jl")



using JSOSolvers, SolverBenchmark, SolverTools, Plots

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


println(" \n\n génération des problemes")
# n_array = [100,500,1000,2000,5000]
n_array = [10,20,30]
# n_array = [1000,2000]
# n_array = [100,200]
problems = create_all_problems(n_array)



println("\n\ndéfinition des solver\n\n")


solver = Dict{Symbol,Function}(
  :trunk => ((prob;kwargs...) -> JSOSolvers.trunk(prob;kwargs...)),
  :trunk_lsr1 => (prob; kwargs...) -> JSOSolvers.trunk(NLPModels.LSR1Model(prob); kwargs...),
  :my_lbfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LBFGS(prob;kwargs...)),
  :my_lsr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.my_LSR1(prob;kwargs...)),
  :p_bfgs => ((prob;kwargs...) -> PartiallySeparableSolvers.PBFGS(prob; kwargs...)),
  :p_sr1 => ((prob;kwargs...) -> PartiallySeparableSolvers.PSR1(prob; kwargs...))
)

const atol = 1.0e-5
const rtol = 1.0e-6
const max_time = 300.0


#= Lancement du benchmark sur les problèmes générés, sur les solvers défini dans la variable solvers =#
println("lancement des benchmmarks")

stats = bmark_solvers(solver, problems; max_time=max_time, max_eval = 5000, atol=atol, rtol=rtol)
keys_solver = keys(stats)

println("affichage du profile des solvers par rapport au problèmes")
# performance_profile(stats, df->df.elapsed_time)
# performance_profile(stats, df->df.iter)

println("affichage des tables")
markdown_table(stdout, stats[:trunk], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad ])
markdown_table(stdout, stats[:trunk_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad ])
markdown_table(stdout, stats[:p_sr1], cols=[  :name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad ])
markdown_table(stdout, stats[:p_bfgs], cols=[ :name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad ])
markdown_table(stdout, stats[:my_lbfgs], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad ])
markdown_table(stdout, stats[:my_lsr1], cols=[:name, :nvar, :elapsed_time, :iter, :dual_feas, :status, :objective, :neval_obj, :neval_grad ])




error("fin")
#= Ecriture des résultats dans un fichier au format markdown=#
println("écriture des résultats markdown")
location_md = string("src/comparaison/results/result_bench_md.txt")
io = open(location_md,"w")
close(io)
io = open(location_md,"w+")
for i in keys_solver
  markdown_table(io, stats[i], cols=[ :name, :nvar, :elapsed_time, :iter, :status, :objective, :neval_obj, :neval_grad ])
end
close(io)


#= Ecriture des résultats dans un fichier au format latex=#
println("écriture des résultats latex")
location_latex = string("src/comparaison/results/result_bench_latex.txt")
io = open(location_latex,"w")
close(io)
io = open(location_latex,"w+")
for i in keys_solver
  latex_table(io, stats[i], cols=[ :name, :nvar, :elapsed_time, :iter, :status, :objective, :neval_obj, :neval_grad ])
end
close(io)










#= non nécessaire actuellement, faisant office de test pour les dictionnaires et la création de fichier de résultat=#
#=
prob_name_file = Dict{Symbol,String}(
        :rosenbrock => "Rosenbrock",
        :chained_wood => "Chained_Wood",
        :chained_Powel => "Chained_Powel",
        :chained_Cragg_Levis => "Chained_Cragg_Levis",
        )

solver_name = Dict{Symbol,String}(
        :psr1 => "P-SR1",
        :trunk => "Trunk",
        :trunk_lsr1 => "Trunk_LSR1",
        :lsr1 => "L-SR1",
        :lbfgs => "L-BFGS",
        :pbfgs => "P-BFGS",
        )

solver_function = Dict{Symbol,String}(
        :psr1 => "P-SR1",
        :trunk => "Trunk",
        :trunk_lsr1 => "Trunk_LSR1",
        :lsr1 => "L-SR1",
        :lbfgs => "L-BFGS",
        :pbfgs => "P-BFGS",
        )

function open_close_all_result_file(prob_name_file, solver_name)
  for (k_prob, name_prob) ∈ prob_name_file
    depo_name = string("src/comparaison/results/",name_prob)
    for (k_solver, name_solver) ∈ solver_name
      file_name = string(name_prob, "_", name_solver,".txt")
      location = string(depo_name,"/", file_name)
      #on ouvre et on ferme le fichier, on le reset
      io = open(location,"w")
      close(io)
    end
  end
end
=#
