# using JuMP, MathOptInterface, LinearAlgebra, SparseArrays
# using Test, BenchmarkTools, ProfileView, InteractiveUtils

using NLPModels, JuMP, MathOptInterface, NLPModelsJuMP


function create_chained_wood_JuMP_Model(n :: Int)
    σ = 10e-5
    m = Model()
    @variable(m, x[1:n])
    @NLobjective(m, Min, sum( 100 * (x[2*j-1]^2 - x[2*j])^2 + (x[2*j-1] - 1)^2 + 90 * (x[2*j+1]^2 - x[2*j+2])^2 + (x[2*j+1] -1)^2 + 10 * (x[2*j] + x[2*j+2] - 2)^2 + (x[2*j] - x[2*j+2])^2 * 0.1  for j in 1:Integer((n-2)/2) )) #rosenbrock function
    evaluator = JuMP.NLPEvaluator(m)
    MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
    obj = MathOptInterface.objective_expr(evaluator)
    vec_var = JuMP.all_variables(m)
    vec_value = create_initial_point_chained_wood(n)
    JuMP.set_start_value.(vec_var, vec_value)
    return (m, evaluator,obj)
end


function create_initial_point_chained_wood(n)
    point_initial = Vector{Float64}(undef, n)
    for i in 1:n
        if i <= 4 && mod(i,2) == 1
            point_initial[i] = -3
        elseif i <= 4 && mod(i,2) == 0
            point_initial[i] = -1
        elseif i > 4 && mod(i,2) == 1
            point_initial[i] = -2
        elseif i > 4 && mod(i,2) == 0
            point_initial[i] = 0
        else
            error("bizarre")
        end
    end
    return point_initial
end

# n_array = [100, 200, 300, 500]
#
# for i in n_array
#     println(" nouveau modèle i = ", i)
#     (m,evaluator,obj) = create_chained_wood_JuMP_Model(i)
#     println("fin de la définition du modèle JuMP")
#     inital_point = create_initial_point_chained_wood(i)
#     println("fin de la définition du point iniitial")
#
#     (cpt2,s) = My_SPS_Model_Module.solver_TR_PSR1!(obj, i, inital_point)
#     nlp = MathOptNLPModel(m)
#     B = LSR1Operator(i, scaling=true) :: LSR1Operator{Float64} #scaling=true
#     (x_f,cpt)  = solver_L_SR1_Ab_NLP(nlp, B, inital_point)
#     @show MathOptInterface.eval_objective(evaluator, inital_point)
#     @show MathOptInterface.eval_objective(evaluator, s.tpl_x[Int(s.index)])
#     @show MathOptInterface.eval_objective(evaluator, x_f)
#     g1 = Vector{Float64}(undef,i)
#     g2 = Vector{Float64}(undef,i)
#     MathOptInterface.eval_objective_gradient(evaluator, s.tpl_x[Int(s.index)],g1)
#     MathOptInterface.eval_objective_gradient(evaluator, x_f, g2)
#     @show norm(g1,2)
#     @show norm(g2,2)
# end
#
#
(m,evaluator,obj_expr) = create_chained_wood_JuMP_Model(8)
# @profview cpt,s = My_SPS_Model_Module.solver_TR_PSR1!(obj_expr, 500, create_initial_point_chained_wood(500))
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
