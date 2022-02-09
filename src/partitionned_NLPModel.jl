using CalculusTreeTools, PartiallySeparableNLPModels
using LinearOperators, NLPModels, JuMP, FastClosures

import NLPModels: increment!, obj,grad!, hprod!, hprod


mutable struct PartionnedNLPModel{T,Y <: Number,S} <: AbstractNLPModel{Y,S}
    meta :: NLPModelMeta{Y,S}
    s_a :: struct_algo{T,Y}
    counters :: Counters
end

function PartionnedNLPModel(nlp :: T) where T <: AbstractNLPModel
    model = nlp.eval.m
    evaluator = JuMP.NLPEvaluator(model)
    MathOptInterface.initialize(evaluator, [:ExprGraph])
    obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
    n = model.moi_backend.model_cache.model.num_variables_created
    x0 = nlp.meta.x0
    return PartionnedNLPModel(obj_Expr, n, x0)
end

function PartionnedNLPModel(nlp :: MathOptNLPModel)
  model = nlp.eval.m
  evaluator = JuMP.NLPEvaluator(model)
  MathOptInterface.initialize(evaluator, [:ExprGraph])
  obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
  n = model.moi_backend.model_cache.model.num_variables_created
  x0 = nlp.meta.x0
  return PartionnedNLPModel(obj_Expr, n, x0)
end

function PartionnedNLPModel(model :: JuMP.Model)
  evaluator = JuMP.NLPEvaluator(model)
  MathOptInterface.initialize(evaluator, [:ExprGraph])
  obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
  n = model.moi_backend.model_cache.model.num_variables_created
  x0 = model.meta.x0
  return PartionnedNLPModel(obj_Expr, n, x0)
end

function PartionnedNLPModel(obj :: T, n :: Int, x0 :: AbstractVector{Y}, t=Y :: DataType) where T where Y <: Number
  meta = NLPModelMeta(n,x0=x0)
  s_a = alloc_struct_algo(obj, n :: Int, t :: DataType )
  return PartionnedNLPModel{CalculusTreeTools.complete_expr_tree,t,Vector{t}}(meta, s_a, Counters())
end

function NLPModels.obj(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector{Y}) where T where Y <: Number
    increment!(nlp, :neval_obj)
    return PartiallySeparableNLPModels.evaluate_SPS(nlp.s_a.sps, x)
end

function NLPModels.grad!(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector{Y}, g :: AbstractVector{Y}) where T where Y <: Number
  increment!(nlp, :neval_grad)
  PartiallySeparableNLPModels.evaluate_SPS_gradient!(nlp.s_a.sps, x, nlp.s_a.tpl_g[1])
  PartiallySeparableNLPModels.build_gradient!(nlp.s_a.sps, nlp.s_a.tpl_g[1], g)
end

function NLPModels.grad(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector{Y}) where T where Y <: Number
    g = similar(x)
    NLPModels.grad!(nlp, x, g)
    return g
end

function NLPModels.hprod!(nlp :: PartionnedNLPModel{T,Y},
                  x :: AbstractVector{Y},
                  v :: AbstractVector{Y},
                  hv :: AbstractVector{Y};
                  obj_weight=1.0,
                  y=Float64[]) where T where Y <: Number
  nlp.counters.neval_hprod += 1
  PartiallySeparableNLPModels.Hv!(hv, nlp.s_a.sps, x, v) 
end                 

function NLPModels.hprod!(nlp :: PartionnedNLPModel{T,Y},                     
                  v :: AbstractVector{Y},
                  hv :: AbstractVector{Y};
                  obj_weight=1.0,
                  y=Float64[]) where T where Y <: Number
  nlp.counters.neval_hprod += 1
  PartiallySeparableNLPModels.Hv!(hv, nlp.s_a.sps, nlp.s_a.sps.x, v)
end 

function NLPModels.hprod(nlp :: PartionnedNLPModel{T,Y},
                  x :: AbstractVector{Y},
                  v :: AbstractVector{Y};
                  obj_weight=1.0,
                  y=Float64[]) where T where Y <: Number
  nlp.counters.neval_hprod += 1
  PartiallySeparableNLPModels.Hv(nlp.s_a.sps, x, v) 
end

function NLPModels.hprod(nlp :: PartionnedNLPModel{T,Y},
                  v :: AbstractVector{Y};
                  obj_weight=1.0,
                  y=Float64[]) where T where Y <: Number
  nlp.counters.neval_hprod += 1
  hv = similar(v) ; PartiallySeparableNLPModels.Hv!(hv, nlp.s_a.sps, nlp.s_a.sps.x, v)
  return hv
end

function NLPModels.hess_op!(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector, Hv :: AbstractVector; obj_weight::Real=one(eltype(x))) where T where Y <: Number
  @lencheck nlp.meta.nvar x Hv
  prod = @closure (res,v) -> NLPModels.hprod!(nlp, v, res)
  return LinearOperator{eltype(x)}(nlp.meta.nvar, nlp.meta.nvar, true, true, prod, prod, prod)
end