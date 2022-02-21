module Mod_partitioned_methods
	using JuMP, MathOptInterface, ModelingToolkit
	using ADNLPModels, NLPModels, NLPModelsJuMP
	using CalculusTreeTools, PartiallySeparableNLPModels
	using ..Mod_TR_CG_part_data
	
	export PBFGS2, PLBFGS

	function get_expr_tree(nlp :: MathOptNLPModel; x0 :: AbstractVector=copy(nlp.meta.x0), kwargs...)
		model = nlp.eval.m
		evaluator = JuMP.NLPEvaluator(model)
		MathOptInterface.initialize(evaluator, [:ExprGraph])
		obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
		ex = CalculusTreeTools.transform_to_expr_tree(obj_Expr)
		n = nlp.meta.nvar
		T = eltype(x0)
		return ex, n, x0, T
	end

	function get_expr_tree(adnlp :: ADNLPModel; x0 :: AbstractVector=copy(adnlp.meta.x0), kwargs...)
		n = adnlp.meta.nvar
		ModelingToolkit.@variables x[1:n]
		fun = adnlp.f(x)
		ex = CalculusTreeTools.transform_to_expr_tree(fun)
		T = eltype(x0)
		return ex, n, x0, T
	end

	function PBFGS2(nlp :: N; kwargs...) where N <: AbstractNLPModel
		(ex, n, x0, T) = get_expr_tree(nlp)
		part_data_pbfgs = build_PartitionedData_TR_PBFGS(ex, n; type=T, x0=x0)
		ges = Generic_algorithm_wrapper(nlp, part_data_pbfgs; kwargs...)
		return ges
	end 

	function PLBFGS(nlp :: N; kwargs...) where N <: AbstractNLPModel
		(ex, n, x0, T) = get_expr_tree(nlp)
		part_data_plbfgs = build_PartitionedData_TR_PLBFGS(ex, n; type=T, x0=x0)
		ges = Generic_algorithm_wrapper(nlp, part_data_plbfgs; kwargs...)
		return ges
end 
end 