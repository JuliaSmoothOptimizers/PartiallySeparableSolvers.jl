function create_initial_point_Rosenbrock(n)
	point_initial = Vector{Float64}(undef, n)
	for i in 1:n
		if mod(i,2) == 1
			point_initial[i] = -1.2
		else mod(i,2) == 0
			point_initial[i] = 1.0
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

function create_Rosenbrock_MathOptModel(n :: Int)
	model_jump , evaluator, obj = create_Rosenbrock_JuMP_Model(n)
	mathoptnlp = MathOptNLPModel(model_jump, name="Rosenbrock "*string(n))
	return mathoptnlp
end 

name_Rosenbrock(n :: Int) = "Rosenbrock_" * string(n)
function Rosenbrock(x)
	n = length(x)
	return sum( 100 * (x[j-1]^2 - x[j])^2 + (x[j-1] - 1)^2  for j in 2:n)
end
Rosenbrock_ADNLPModel(n :: Int) = ADNLPModel(eval(Rosenbrock), create_initial_point_Rosenbrock(n), name=name_Rosenbrock(n))
