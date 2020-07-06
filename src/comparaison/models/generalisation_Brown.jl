using NLPModels, JuMP, MathOptInterface, NLPModelsJuMP


"""
    generalisation_Brown(n)
Create a Brown function model using the package JuMP. It return (m,evaluator,obj), evaluator is usefull to use MOI function
and obj is the Expr that describe the objective function of m.
"""
function create_generalisation_Brown(n :: Int)
    m = Model()
    @variable(m, x[1:n])
    @NLobjective(m, Min, sum( (x[j-1]-3)^2 + (x[j-1]-x[j])^2 + exp(20 * (x[j-1] - x[j]))  for j in 2:n) )  #chained powel
    evaluator = JuMP.NLPEvaluator(m)
    MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
    obj = MathOptInterface.objective_expr(evaluator)
    vec_var = JuMP.all_variables(m)
    vec_value = create_initial_point_chained_cragg_levy(n)
    JuMP.set_start_value.(vec_var, vec_value)
    return (m, evaluator,obj)
end

"""
    create_initial_point_chained_cragg_levy(n)
Create the initial point of Cragg Levy chained function according to the article of Luksan & Vlcek. the initial point of size n
"""
function create_initial_point_brown(n)
    point_initial = Vector{Float64}(undef, n)
    for i in 1:n
        point_initial[i] = 1.0
    end
    return point_initial
end

name_model_brown(n :: Int) = "Brown " * string(n)

function generalisation_brown(x)
  n = length(x)
  return sum( (x[i-1] - 3 )^2 + (x[i-1] - x[i])^2 + exp(20 * (x[i-1] - x[i]))  for i=2:n)
end



create_generalisation_brown_ADNLPModel(n :: Int) = RADNLPModel(eval(generalisation_brown), create_initial_point_brown(n), name=name_model_brown(n))
function create_generalisation_brown_JuMPModel(n :: Int)
    (m_brown, evaluator, obj) = create_generalisation_Brown(n)
    return MathOptNLPModel(m_brown, name=name_model_brown(n))
end
