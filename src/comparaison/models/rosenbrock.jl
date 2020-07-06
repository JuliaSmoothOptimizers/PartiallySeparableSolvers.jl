using NLPModels, JuMP, MathOptInterface, NLPModelsJuMP
using Printf


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


name_model_rosenbrock(n :: Int) = "Rosenbrock " * string(n)

function rosenbrock(x)
  n = length(x)
  return sum( 100 * (x[i-1]^2 - x[i])^2 + (x[i-1] -1 )^2  for i=2:n)
end



create_rosenbrock_ADNLPModel(n :: Int) = RADNLPModel(eval(rosenbrock), create_initial_point_Rosenbrock(n), name=name_model_rosenbrock(n))
function create_rosenbrock_JuMPModel(n :: Int)
    (m_ros, evaluator, obj) = create_Rosenbrock_JuMP_Model(n)
    return MathOptNLPModel(m_ros, name=name_model_rosenbrock(n))
end
