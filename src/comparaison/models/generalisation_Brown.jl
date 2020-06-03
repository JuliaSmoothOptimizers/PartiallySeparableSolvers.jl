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
function create_initial_point_chained_cragg_levy(n)
    point_initial = Vector{Float64}(undef, n)
    for i in 1:n
        point_initial[i] = 1.0
    end
    return point_initial
end



 # (m, evaluator,obj) = create_generalisation_Brown(8)
