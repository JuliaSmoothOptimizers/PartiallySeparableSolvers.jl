using JuMP, MathOptInterface, NLPModelsJuMP, LinearAlgebra, SparseArrays, Test
using CalculusTreeTools, PartiallySeparableNLPModel, PartiallySeparableSolvers
using JSOSolvers

using NLPModels

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
Struct_PS = PartiallySeparableNLPModel.deduct_partially_separable_structure(obj, n)

JuMP_mod = MathOptNLPModel(m_ros, name="Ros "*string(n))

SPS_nlp = PartiallySeparableSolvers.PartionnedNLPModel(JuMP_mod)
x = (x -> 2*x).(ones(SPS_nlp.meta.nvar))
v = ones(SPS_nlp.meta.nvar)


@testset "obj/grad/hprod" begin


    obj_sps = NLPModels.obj(SPS_nlp, x)
    obj_jump = NLPModels.obj(JuMP_mod, x)
    @test obj_jump == obj_sps


    grad_sps = similar(x)
    grad_jump = similar(x)
    NLPModels.grad!(SPS_nlp, x, grad_sps)
    NLPModels.grad!(JuMP_mod, x, grad_jump)
    @test grad_jump == grad_sps


    hv_sps = similar(x)
    hv_jump = similar(x)
    NLPModels.hprod!(SPS_nlp, x, v, hv_sps)
    NLPModels.hprod!(JuMP_mod, x, v, hv_jump)
    @test hv_sps == hv_jump
end


temp_sps = similar(x)
temp_jump = similar(x)
H_sps = hess_op!(SPS_nlp, x, temp_sps)
H_jump = hess_op!(JuMP_mod, x, temp_jump)

H_sps * v
@test H_sps * v == H_jump * v
Juno.@run H_sps * v
