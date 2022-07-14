module Mod_partitioned_methods
using JuMP, MathOptInterface, ModelingToolkit
using ADNLPModels, NLPModels, NLPModelsJuMP
using ExpressionTreeForge, PartiallySeparableNLPModels
using ..Mod_TR_CG_part_data

export get_expr_tree
export PUS

"""
    expr_tree, n, x0 = get_expr_tree(nlp::MathOptNLPModel; x0::Vector{T} = copy(nlp.meta.x0), kwargs...) where {T <: Number}
    expr_tree, n, x0 = get_expr_tree(adnlp::ADNLPModel; x0::Vector{T} = copy(adnlp.meta.x0), kwargs...) where {T <: Number}

Return the `expr_tree`, the size `n` and the initial point `x0` from either a `MathOptNLPModel` or a `ADNLPModel`.
"""
function get_expr_tree(
  nlp::MathOptNLPModel;
  x0::Vector{T} = copy(nlp.meta.x0),
  kwargs...,
) where {T <: Number}
  model = nlp.eval.m
  evaluator = JuMP.NLPEvaluator(model)
  MathOptInterface.initialize(evaluator, [:ExprGraph])
  obj_Expr = MathOptInterface.objective_expr(evaluator)::Expr
  expr_tree =
    ExpressionTreeForge.transform_to_expr_tree(obj_Expr)::ExpressionTreeForge.Type_expr_tree
  n = nlp.meta.nvar
  return expr_tree, n, x0
end

function get_expr_tree(
  adnlp::ADNLPModel;
  x0::Vector{T} = copy(adnlp.meta.x0),
  kwargs...,
) where {T <: Number}
  n = adnlp.meta.nvar
  ModelingToolkit.@variables x[1:n]
  fun = adnlp.f(x)
  expr_tree = ExpressionTreeForge.transform_to_expr_tree(fun)::ExpressionTreeForge.Type_expr_tree
  return expr_tree, n, x0
end

"""
    ges = PUS(nlp::AbstractNLPModel; name = :plse, kwargs...)

`PUS` (partitioned update solver) return a `ges::GenericExecutionStats` from a partitioned quasi-Newton trust-region method.
It can perform several variants, which can be select with the optional argument `name`.
You will perform a:

* PBFGS method with `name=:pbfgs`;
* PSR1 method with `name=:psr1`;
* PSE method with `name=:pse`;
* PLBFGS method with `name=:plbfgs`;
* PLSR1 method with `name=:plsr1`;
* PLSE method with `name=:plse`, by default.
"""
function PUS(nlp::AbstractNLPModel; name = :plse, kwargs...)
  (ex, n, x0) = get_expr_tree(nlp)
  part_data_pqn = build_PartitionedDataTRPQN(ex, n; name = name, x0 = x0)
  ges = generic_algorithm_wrapper(nlp, part_data_pqn; kwargs...)
  return ges
end

end
