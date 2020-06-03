# using JuMP, MathOptInterface, LinearAlgebra, SparseArrays
# using Test, BenchmarkTools, ProfileView, InteractiveUtils

using NLPModels, JuMP, MathOptInterface, NLPModelsJuMP
using Printf

# include("../ordered_include.jl")








function create_Rosenbrock_JuMP_Model(n :: Int)
    m = Model()
    @variable(m, x[1:n])
    @NLobjective(m, Min, sum( 100 * (x[j-1]^2 - x[j])^2 + (x[j-1] - 1)^2  for j in 2:n)) #rosenbrock function
    evaluator = JuMP.NLPEvaluator(m)
    MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
    obj = MathOptInterface.objective_expr(evaluator)
    vec_var = JuMP.all_variables(m)
    vec_value = create_initial_point_Rosenbrock(n)
    JuMP.set_start_value.(vec_var, vec_value)
    return (m, evaluator,obj)
end


function create_initial_point_Rosenbrock(n)
    point_initial = Vector{Float64}(undef, n)
    for i in 1:n
        if mod(i,2) == 1
            point_initial[i] = -1.2
        elseif mod(i,2) == 0
            point_initial[i] = 1.0
        else
            error("bizarre")
        end
    end
    return point_initial
end

# (m,evaluator,obj) = create_Rosenbrock_JuMP_Model(50)
# nlp = MathOptNLPModel(m)
# @show nlp.meta.x0


# n_array = [100, 500, 1000, 2000, 3000, 5000]
# n_array = [100,500,1000, 2000, 3000, 5000, 10000, 20000]
#
#
# for i in n_array
#     println(" \n\n nouveau modèle à ", i, " variables")
#     (m,evaluator,obj) = create_Rosenbrock_JuMP_Model(i)
#     println("fin de la définition du modèle JuMP")
#     initial_point = create_initial_point_Rosenbrock(i)
#     println("fin de la définition du point iniitial")
#
#     t1 = @elapsed ((cpt,s) = My_SPS_Model_Module.solver_TR_PSR1!(obj, i, initial_point))
#
#     nlp = MathOptNLPModel(m)
#     B = LSR1Operator(i, scaling=true) :: LSR1Operator{Float64} #scaling=true
#     t2 = @elapsed (x_f,cpt2)  = solver_L_SR1_Ab_NLP(nlp, B, initial_point)
#
#     println("\nla valeur de la fonction au point initial est ", MathOptInterface.eval_objective(evaluator, initial_point))
#     println("Pour la méthode PSR1 on a fait ", cpt, " itérations pour trouver une valeur de fonction objectif de ", MathOptInterface.eval_objective(evaluator, s.tpl_x[Int(s.index)]), " en ",  t1, " secondes et donc une moyenne de ", t1/cpt, "seconde par itérations")
#     println("Pour la méthode LSR1 on a fait ", cpt2, " itérations pour trouver une valeur de fonction objectif de ", MathOptInterface.eval_objective(evaluator, x_f), " en ",  t2, " secondes et donc une moyenne de ", t2/cpt2, "seconde par itérations")
#     println("le rapport iteration PSR1/iteration LSR1 est ", (t1/cpt)/(t2/cpt2))
# end
#
# println("fin de la boucle")
# println("résolution PSR1")
#
# println("résolution LSR1 JuMP")
#
#

# @show MathOptInterface.eval_objective(evaluator, point_initial)
# @show MathOptInterface.eval_objective(evaluator, x_f)
# @show MathOptInterface.eval_objective(evaluator, s.tpl_x[Int(s.index)])
# @show MathOptInterface.eval_objective_gradient(evaluator, s.tpl_x[Int(s.index)])




# i = 1000
# (m,evaluator,obj) = create_Rosenbrock_JuMP_Model(i)
# println("fin de la définition du modèle JuMP")
# initial_point = create_initial_point_Rosenbrock(i)
# println("fin de la définition du point iniitial")
# #
# t1 = @elapsed ((cpt,s) = My_SPS_Model_Module.solver_TR_PSR1!(obj, i, initial_point))

# val, t, bytes, gctime, memallocs = @timed  My_SPS_Model_Module.solver_TR_PSR1!(obj, i, initial_point)
# error("fin anticipé")
#
