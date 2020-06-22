
using JuMP, MathOptInterface, NLPModelsJuMP, NLPModels
using LinearAlgebra, SparseArrays, Test
using CalculusTreeTools, PartiallySeparableNLPModel, PartiallySeparableSolvers

function create_chained_wood_JuMP_Model(n :: Int)
    σ = 10e-5
    m = Model()
    @variable(m, x[1:n])
    @NLobjective(m, Min, sum( 100 * (x[2*j-1]^2 - x[2*j])^2 + (x[2*j-1] - 1)^2 + 90 * (x[2*j+1]^2 - x[2*j+2])^2 + (x[2*j+1] -1)^2 + 10 * (x[2*j] + x[2*j+2] - 2)^2 + (x[2*j] - x[2*j+2])^2 * 0.1  for j in 1:Integer((n-2)/2) )) #rosenbrock function
    evaluator = JuMP.NLPEvaluator(m)
    MathOptInterface.initialize(evaluator, [:ExprGraph, :Hess])
    obj = MathOptInterface.objective_expr(evaluator)
    vec_var = JuMP.all_variables(m)
    vec_value = create_initial_point_chained_wood(n)
    JuMP.set_start_value.(vec_var, vec_value)
    return (m, evaluator,obj)
end


function create_initial_point_chained_wood(n)
    point_initial = Vector{Float64}(undef, n)
    for i in 1:n
        if i <= 4 && mod(i,2) == 1
            point_initial[i] = -3
        elseif i <= 4 && mod(i,2) == 0
            point_initial[i] = -1
        elseif i > 4 && mod(i,2) == 1
            point_initial[i] = -2
        elseif i > 4 && mod(i,2) == 0
            point_initial[i] = 0
        else
            error("bizarre")
        end
    end
    return point_initial
end

n = 100
(m_ros,evaluator,obj) = create_chained_wood_JuMP_Model(n)
obj_expr_tree = CalculusTreeTools.transform_to_expr_tree(obj)

complete_expr_tree = CalculusTreeTools.create_complete_tree(obj_expr_tree)
Struct_PS = PartiallySeparableNLPModel.deduct_partially_separable_structure(complete_expr_tree, n)

JuMP_mod = MathOptNLPModel(m_ros, name="Chained Powel "*string(n))

#vérification que certaines fonction objectifs sont convexes et d'autre inconnus
# @show Struct_PS.structure[1]
# @show Struct_PS.structure[2]
# @show Struct_PS.structure[3]
# @show Struct_PS.structure[4]
# @show Struct_PS.structure[5]
# @show Struct_PS.structure[6]

@testset " test convexity detection chainedwood" begin
    for i in 1:length(Struct_PS.structure)
        if (i%6 == 1) || (i%6 == 3)
            @test CalculusTreeTools.is_unknown(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[i]))
        else
            @test CalculusTreeTools.is_convex(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[i]))
        end
    end

    # @test CalculusTreeTools.is_unknown(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[1]))
    # @test CalculusTreeTools.is_convex(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[2]))
    # @test CalculusTreeTools.is_unknown(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[3]))
    # @test CalculusTreeTools.is_convex(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[4]))
    # @test CalculusTreeTools.is_convex(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[5]))
    # @test CalculusTreeTools.is_convex(PartiallySeparableNLPModel.get_convexity_status(Struct_PS.structure[6]))
end
