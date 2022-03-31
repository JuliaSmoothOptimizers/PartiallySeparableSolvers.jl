module Mod_partitioned_methods
	using JuMP, MathOptInterface, ModelingToolkit
	using ADNLPModels, NLPModels, NLPModelsJuMP
	using CalculusTreeTools, PartiallySeparableNLPModels
	using ..Mod_TR_CG_part_data
	
	export get_expr_tree
	export PBFGS2, PLBFGS, PUS

	function get_expr_tree(nlp :: MathOptNLPModel; x0 :: Vector{T}=copy(nlp.meta.x0), kwargs...) where T <: Number
		model = nlp.eval.m
		evaluator = JuMP.NLPEvaluator(model)
		MathOptInterface.initialize(evaluator, [:ExprGraph])
		obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
		ex = CalculusTreeTools.transform_to_expr_tree(obj_Expr) :: CalculusTreeTools.t_expr_tree
		n = nlp.meta.nvar		
		return ex, n, x0
	end

	function get_expr_tree(adnlp :: ADNLPModel; x0 :: Vector{T}=copy(adnlp.meta.x0), kwargs...) where T <: Number
		n = adnlp.meta.nvar
		ModelingToolkit.@variables x[1:n]
		fun = adnlp.f(x)
		ex = CalculusTreeTools.transform_to_expr_tree(fun) :: CalculusTreeTools.t_expr_tree		
		return ex, n, x0
	end

	function PBFGS2(nlp :: N; kwargs...) where N <: AbstractNLPModel
		(ex, n, x0) = get_expr_tree(nlp) 
		part_data_pbfgs = build_PartitionedData_TR_PBFGS(ex, n; x0=x0) 
		ges = generic_algorithm_wrapper(nlp, part_data_pbfgs; kwargs...)
		return ges
	end 

	function PLBFGS(nlp :: N; kwargs...) where N <: AbstractNLPModel
		(ex, n, x0) = get_expr_tree(nlp)
		part_data_plbfgs = build_PartitionedData_TR_PLBFGS(ex, n; x0=x0)
		ges = generic_algorithm_wrapper(nlp, part_data_plbfgs; kwargs...)
		return ges
	end 

	function PUS(nlp :: N; name=:plse, kwargs...) where N <: AbstractNLPModel
		(ex, n, x0) = get_expr_tree(nlp)
		part_data_pqn = build_PartitionedData_TR_PQN(ex, n; name=name, x0=x0)
		ges = generic_algorithm_wrapper(nlp, part_data_pqn; kwargs...)
		return ges
	end 
	
end 