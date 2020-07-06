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



name_model_cragg_levy(n :: Int) = "CraggLevy " * string(n)

function cragglevy(x)
  n = length(x)
  (n-2) % 2 == 0 || error("number of variables minus 2 must be even")
  return sum( exp(x[2*i-1])^4 + 100 * (x[2*i]- x[2*i+1])^6 + (tan(x[2*i+1]) - x[2*i+2])^4 + x[2*i-1]^8 + (x[2*i+2] -1 )^2 for i=1:div(n-2,2))
end


create_cragg_levy_ADNLPModel(n :: Int) = RADNLPModel(eval(cragglevy), create_initial_point_chained_cragg_levy(n), name=name_model_cragg_levy(n))
function create_cragg_levy_JuMPModel(n :: Int)
  (m_cragg, evaluator,obj) = create_chained_cragg_levy_JuMP_Model(n)
  return MathOptNLPModel(m_cragg, name=name_model_cragg_levy(n))
end
