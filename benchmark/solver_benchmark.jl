using BenchmarkTools
using JSOSolvers, SolverTools
using NLPModelsJuMP, LinearOperators, NLPModels
using JuMP, MathOptInterface
using PartiallySeparableSolvers


const SUITE = BenchmarkGroup()


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


n = [100,200,500,1000,2000,5000]
# n = [10,20,30]
# n = [500,1000,2000,4000]
problems = create_Rosenbrock_JuMP_Model.(n)
nlp_problems = MathOptNLPModel.([p[1] for p in problems])



SUITE["Trunk"] = BenchmarkGroup()
SUITE["Trunk_LSR1"] = BenchmarkGroup()
SUITE["L-BFGS"] = BenchmarkGroup()
SUITE["L-SR1"] = BenchmarkGroup()
SUITE["P-SR1"] = BenchmarkGroup()
SUITE["P-BFGS"] = BenchmarkGroup()
SUITE["P-BS"] = BenchmarkGroup()

for i in 1:length(problems)

  atol = 1.0e-5
  rtol = 1.0e-6
  max_time = 300.0
  max_eval = 5000
  (m_ros, evaluator, obj_ros) = problems[i]
  n = m_ros.moi_backend.model_cache.model.num_variables_created
  prob = nlp_problems[i]
  LSR1_prob = NLPModels.LSR1Model(prob)

  SUITE["ros $n var"] = BenchmarkGroup()

  SUITE["ros $n var"]["Trunk"] = @benchmarkable $(JSOSolvers.trunk)($prob)
  SUITE["ros $n var"]["Trunk_LSR1"] = @benchmarkable $(JSOSolvers.trunk)($LSR1_prob)
  SUITE["ros $n var"]["L-BFGS"] = @benchmarkable $(PartiallySeparableSolvers.my_LBFGS)($prob)
  SUITE["ros $n var"]["L-SR1"] = @benchmarkable $(PartiallySeparableSolvers.my_LSR1)($prob)
  SUITE["ros $n var"]["P-BFGS"] = @benchmarkable $(PartiallySeparableSolvers.PBFGS)($prob)
  SUITE["ros $n var"]["P-SR1"] = @benchmarkable $(PartiallySeparableSolvers.PSR1)($prob)
  SUITE["ros $n var"]["P-BS"] = @benchmarkable $(PartiallySeparableSolvers.PBS)($prob)

end
