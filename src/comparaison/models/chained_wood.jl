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

name_model_chained_wood(n :: Int) = "ChainedWood " * string(n)

function chainwoo(x)
  n = length(x)
  n % 4 == 0 || error("number of variables must be a multiple of 4")
  return 1.0 + sum(100 * (x[2*i]   - x[2*i-1]^2)^2 + (1 - x[2*i-1])^2 +
              90 * (x[2*i+2] - x[2*i+1]^2)^2 + (1 - x[2*i+1])^2 +
              10 * (x[2*i] + x[2*i+2] - 2)^2 + 0.1 * (x[2*i] - x[2*i+2])^2 for i=1:div(n,2)-1)
end


create_chained_wood_ADNLPModel(n :: Int) = ADNLPModel(eval(chainwoo), create_initial_point_chained_wood(n), name=name_model_chained_wood(n))
function create_chained_wood_JuMPModel(n :: Int)
  (m_chained, evaluator,obj) = create_chained_wood_JuMP_Model(n)
  return MathOptNLPModel(m_chained, name=name_model_chained_wood(n))
end
