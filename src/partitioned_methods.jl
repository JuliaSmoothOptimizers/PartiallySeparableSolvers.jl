module ModPartitionedMethods
using JuMP, MathOptInterface, ModelingToolkit
using ADNLPModels, NLPModels, NLPModelsJuMP
using ExpressionTreeForge, PartiallySeparableNLPModels
using ..ModTrustRegionPartitionedData

export get_expr_tree
export PTRUNK

"""
    expr_tree, n, x0 = get_expr_tree(nlp::MathOptNLPModel; x0::Vector{T} = copy(nlp.meta.x0), kwargs...) where {T <: Number}
    expr_tree, n, x0 = get_expr_tree(adnlp::ADNLPModel; x0::Vector{T} = copy(adnlp.meta.x0), kwargs...) where {T <: Number}

Return the `expr_tree`, the size `n` and the initial point `x0` from either a `MathOptNLPModel` or a `ADNLPModel`.
"""
function get_expr_tree(
  nlp::MathOptNLPModel;
  x0::AbstractVector{T} = copy(nlp.meta.x0),
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
  x0::AbstractVector{T} = copy(adnlp.meta.x0),
  kwargs...,
) where {T <: Number}
  n = adnlp.meta.nvar
  ModelingToolkit.@variables x[1:n]
  fun = adnlp.f(x)
  expr_tree = ExpressionTreeForge.transform_to_expr_tree(fun)::ExpressionTreeForge.Type_expr_tree
  return expr_tree, n, x0
end

"""
    stats = PTRUNK(nlp::AbstractNLPModel; name = :plse, kwargs...)

`PTRUNK` (partitioned update solver) return a `stats::GenericExecutionStats` from a partitioned quasi-Newton trust-region method.
Several variants are available via the optional argument `name`.
Each variant approximate the element-Hessians differently, some use dense matrices while others use linear-operators.
Variants with dense element-Hessian approximations:

* PBFGS with `name=:pbfgs`, every element-Hessian approximation is updated with BFGS;
* PSR1 with `name=:psr1`, every element-Hessian approximation is updated with SR1;
* PSE with `name=:pse`, every element-Hessian approximation is updated with BFGS or SR1 if the curvature condition doesn't hold.

Variants with linear-operator element-Hessian approximations:

* PLBFGS with `name=:plbfgs`, every element-Hessian approximations is a LBFGS operator;
* PLSR1 with `name=:plsr1`, every element-Hessian approximations is a LSR1 operator;
* PLSE with `name=:plse`, by default, every element-Hessian approximations is a LBFGS operator as long as the curvature condition holds, otherwise it becomes a LSR1 operator.
"""
function PTRUNK(nlp::AbstractNLPModel; name = :plse, kwargs...)
  (ex, n, x0) = get_expr_tree(nlp)
  part_data_pqn = build_PartitionedDataTRPQN(ex, n; name = name, x0 = x0)
  stats = partitionedTrunk(nlp, part_data_pqn; kwargs...)
  return stats
end

end
