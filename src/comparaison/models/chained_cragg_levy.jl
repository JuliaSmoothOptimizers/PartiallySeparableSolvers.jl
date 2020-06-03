using NLPModels, JuMP, MathOptInterface, NLPModelsJuMP


"""
    create_chained_cragg_levy_JuMP_Model(n)
Create a Chained Cragg Levy Model m using the package JuMP. It return (m,evaluator,obj), evaluator is usefull to use MOI function
and obj is the Expr that describe the objective function of m.
"""
function create_chained_cragg_levy_JuMP_Model(n :: Int)
    m = Model()
    @variable(m, x[1:n])
    @NLobjective(m, Min, sum( exp(x[2*j-1])^4 + 100*(x[2*j] - x[2*j+1])^6 + (tan(x[2*j+1]) - x[2*j+2])^4 + x[2*j-1]^8 + (x[2*j+2]-1)^2 for j in 1:Integer((n-2)/2) ) ) #chained powel
    evaluator = JuMP.NLPEvaluator(m)
    MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
    obj = MathOptInterface.objective_expr(evaluator)
    vec_var = JuMP.all_variables(m)
    vec_value = create_initial_point_chained_cragg_levy(n)
    JuMP.set_start_value.(vec_var, vec_value)
    # JuMP.set_name(m, "Cragg Levy "*string(n))
    # m.moi_backend.model_cache.model.name = "Cragg Levy "*string(n)
    return (m, evaluator,obj)
end

"""
    create_initial_point_chained_cragg_levy(n)
Create the initial point of Cragg Levy chained function according to the article of Luksan & Vlcek. the initial point of size n
"""
function create_initial_point_chained_cragg_levy(n)
    point_initial = Vector{Float64}(undef, n)
    point_initial[1] = 1.0
    for i in 2:n
        point_initial[i] = 2.0
    end
    return point_initial
end



#  (m, evaluator,obj) = create_chained_cragg_levy_JuMP_Model(8)
# @show m.moi_backend.model_cache.model.name
