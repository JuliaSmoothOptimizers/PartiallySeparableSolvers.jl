using CalculusTreeTools, PartiallySeparableNLPModel

using LinearOperators
import NLPModels: increment!, obj,grad!, hprod!, hprod
using NLPModels

using FastClosures

    mutable struct PartionnedNLPModel{T,Y <: Number} <: AbstractNLPModel
        meta :: NLPModelMeta
        s_a :: struct_algo{T, Y}
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

    function PartionnedNLPModel(obj :: T, n :: Int, x0 :: AbstractVector{Y}, t=Y :: DataType) where T where Y <: Number
        meta = NLPModelMeta(n,x0=x0)
        temp_tree = CalculusTreeTools.transform_to_expr_tree(obj)
        work_obj = CalculusTreeTools.create_complete_tree(temp_tree)

        s_a = alloc_struct_algo(work_obj, n :: Int, t :: DataType )

        return PartionnedNLPModel{CalculusTreeTools.complete_expr_tree,t}(meta, s_a, Counters())
    end


    function NLPModels.obj(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector{Y}) where T where Y <: Number
        increment!(nlp, :neval_obj)
        return PartiallySeparableNLPModel.evaluate_SPS(nlp.s_a.sps, x)
    end

    function NLPModels.grad!(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector{Y}, g :: AbstractVector{Y}) where T where Y <: Number
      increment!(nlp, :neval_grad)
      PartiallySeparableNLPModel.evaluate_SPS_gradient!(nlp.s_a.sps, x, nlp.s_a.tpl_g[1])
      PartiallySeparableNLPModel.build_gradient!(nlp.s_a.sps, nlp.s_a.tpl_g[1], g)
    end

    function NLPModels.grad(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector{Y}) where T where Y <: Number
        g = similar(x)
        NLPModels.grad!(nlp, x, g)
        return g
    end

    # function NLPModels.hess(nlp :: SPS_Model{T,Y}, x :: AbstractVector) where T
    #   increment!(nlp, :neval_hess)
    #   return Array(PartiallySeparableStructure.evaluate_hessian(nlp.storage, x ) )
    # end


    NLPModels.hprod!(nlp :: PartionnedNLPModel{T,Y},
                     x :: AbstractVector{Y},
                     v :: AbstractVector{Y},
                     hv :: AbstractVector{Y};
                     obj_weight=1.0,
                     y=Float64[]) where T where Y <: Number = begin nlp.counters.neval_hprod += 1 ; PartiallySeparableNLPModel.Hv!(hv, nlp.s_a.sps, x, v) end

     NLPModels.hprod(nlp :: PartionnedNLPModel{T,Y},
                      x :: AbstractVector{Y},
                      v :: AbstractVector{Y};
                      obj_weight=1.0,
                      y=Float64[]) where T where Y <: Number = begin nlp.counters.neval_hprod += 1 ; PartiallySeparableNLPModel.Hv(nlp.s_a.sps, x, v) end


     NLPModels.hprod(nlp :: PartionnedNLPModel{T,Y},
                      v :: AbstractVector{Y};
                      obj_weight=1.0,
                      y=Float64[]) where T where Y <: Number = begin nlp.counters.neval_hprod += 1; hv = similar(v) ; PartiallySeparableNLPModel.Hv!(hv, nlp.s_a.sps, nlp.s_a.sps.x, v); return hv end



    function NLPModels.hess_op!(nlp :: PartionnedNLPModel{T,Y}, x :: AbstractVector, Hv :: AbstractVector; obj_weight::Real=one(eltype(x))) where T where Y <: Number
        @lencheck nlp.meta.nvar x Hv
        prod = @closure v -> NLPModels.hprod(nlp, x, v; obj_weight=obj_weight)
        return LinearOperator{eltype(x)}(nlp.meta.nvar, nlp.meta.nvar, true, true, prod, prod, prod)
    end
