module ModPartitionedMethods

using ADNLPModels, NLPModels, NLPModelsJuMP
using ExpressionTreeForge 
# PartiallySeparableNLPModels
using ..ModTrustRegionPartitionedData
using ..Mod_PQN

export PTRUNK

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
  x0 = nlp.meta.x0
  n = nlp.meta.nvar
  ex = ExpressionTreeForge.get_expression_tree(nlp)
  part_data_pqn = build_PartitionedDataTRPQN(ex, n; name = name, x0 = x0)
  stats = partitionedTrunk(nlp, part_data_pqn; kwargs...)
  return stats
end

end
