using JuMP, MathOptInterface, NLPModelsJuMP
using CalculusTreeTools, PartiallySeparableNLPModel

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


n = 100
(m_ros,evaluator,obj) = create_Rosenbrock_JuMP_Model(n)
obj_expr_tree = CalculusTreeTools.transform_to_expr_tree(obj)

JuMP_mod = MathOptNLPModel(m_ros, name="Ros "*string(n))

x = create_initial_point_Rosenbrock(n)

# PartiallySeparableSolvers.solver_TR_PSR1!(obj, n, x)
# PartiallySeparableSolvers.solver_TR_PSR1!(obj_expr_tree, n, x)
# PartiallySeparableSolvers.solver_TR_PBFGS!(obj, n, x)
# PartiallySeparableSolvers.solver_TR_PBFGS!(obj_expr_tree, n, x)

ges1 = solver_TR_PBFGS!(JuMP_mod)
ges2 = solver_TR_PSR1!(JuMP_mod)

@test ges1.objective < 1e-5 && ges2.objective < 4
# MathOptInterface.eval_objective(evaluator, ges1.solution)
# MathOptInterface.eval_objective(evaluator, ges2.solution)
# MathOptInterface.eval_objective_gradient()
