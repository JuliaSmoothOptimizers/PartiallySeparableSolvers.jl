using JuMP, MathOptInterface, NLPModelsJuMP, LinearAlgebra, SparseArrays
using CalculusTreeTools, PartiallySeparableNLPModel

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

x = create_initial_point_Rosenbrock(n)
σ = 1e-5

@testset  "vérification des mises à jour SR1/BFGS " begin
    x = ones(n)
    x_1 = (tmp -> 2*tmp).(x)

    f_approx = ( elm_fun :: PartiallySeparableNLPModel.element_function{} -> PartiallySeparableNLPModel.element_hessian{Float64}( Array{Float64,2}(zeros(Float64, length(elm_fun.used_variable), length(elm_fun.used_variable)) ) ) )
    f = (y :: PartiallySeparableNLPModel.element_function{} -> PartiallySeparableNLPModel.element_gradient{typeof(x[1])}(Vector{typeof(x[1])}(undef, length(y.used_variable) )) )

    exact_Hessian = PartiallySeparableNLPModel.Hess_matrix{Float64}(f_approx.(Struct_PS.structure))
    approx_hessian_SR1 = PartiallySeparableNLPModel.Hess_matrix{Float64}(f_approx.(Struct_PS.structure))
    approx_hessian_BFGS = PartiallySeparableNLPModel.Hess_matrix{Float64}(f_approx.(Struct_PS.structure))

    grad_x = PartiallySeparableNLPModel.grad_vector{typeof(x[1])}( f.(Struct_PS.structure) )
    grad_x_1 = PartiallySeparableNLPModel.grad_vector{typeof(x[1])}( f.(Struct_PS.structure) )
    grad_diff = PartiallySeparableNLPModel.grad_vector{typeof(x[1])}( f.(Struct_PS.structure) )


    s = x_1 - x
    PartiallySeparableNLPModel.struct_hessian!(Struct_PS, x, exact_Hessian)
    PartiallySeparableNLPModel.evaluate_SPS_gradient!(Struct_PS, x, grad_x)
    PartiallySeparableNLPModel.evaluate_SPS_gradient!(Struct_PS, x_1, grad_x_1)
    PartiallySeparableNLPModel.minus_grad_vec!(grad_x_1, grad_x, grad_diff)

    #calcul de l'approximation
    PartiallySeparableSolvers.update_SPS_SR1!(Struct_PS, exact_Hessian, approx_hessian_SR1, grad_diff, s)
    PartiallySeparableSolvers.update_SPS_SR1!(Struct_PS, exact_Hessian, approx_hessian_BFGS, grad_diff, s)
    # mettre les informations sous des formats comparable (Vector)
    dif_gradient = PartiallySeparableNLPModel.build_gradient(Struct_PS, grad_diff)
    Bs_SR1 = PartiallySeparableNLPModel.product_matrix_sps(Struct_PS, approx_hessian_SR1, s)
    Bs_BFGS = PartiallySeparableNLPModel.product_matrix_sps(Struct_PS, approx_hessian_BFGS, s)

    check_is_pos_def = (elmt_hess -> isposdef(elmt_hess.elmt_hess) )
    my_and = ( (x,y) -> x && y )
    @test mapreduce( check_is_pos_def, my_and, approx_hessian_BFGS.arr)
    #ce test à voir
    sp_approx_BFGS = PartiallySeparableNLPModel.construct_Sparse_Hessian(Struct_PS, approx_hessian_BFGS)
    @test isposdef(sp_approx_BFGS)

    #test
    @test norm(Bs_BFGS - dif_gradient, 2) < σ
    @test norm(Bs_SR1 - dif_gradient, 2) < σ
end


# PartiallySeparableSolvers.solver_TR_PSR1!(obj, n, x)
# PartiallySeparableSolvers.solver_TR_PSR1!(obj_expr_tree, n, x)
# PartiallySeparableSolvers.solver_TR_PBFGS!(obj, n, x)
# PartiallySeparableSolvers.solver_TR_PBFGS!(obj_expr_tree, n, x)

ges1 = solver_TR_PBFGS!(JuMP_mod)
ges2 = solver_TR_PSR1!(JuMP_mod)
ges4 = my_LSR1(JuMP_mod)
ges3 = my_LBFGS(JuMP_mod)

ges = [ges1,ges2,ges3,ges4]
MOI_gradient = Vector{typeof(x[1])}(undef,n)
ges_grad = []
for i in ges
    MathOptInterface.eval_objective_gradient(evaluator, MOI_gradient, i.solution)
    push!(ges_grad, copy(MOI_gradient))
end
ges_nrm_grad = ( x -> norm(x)).(ges_grad)

grad_init = Vector{typeof(x[1])}(undef,n)
MathOptInterface.eval_objective_gradient(evaluator, grad_init, create_initial_point_Rosenbrock(n))
check_nrm = (nrm_grad -> nrm_grad < 1e-6 * norm(grad_init) )


@test ges1.objective < 1e-5 && ges2.objective < 4 && ges3.objective < 1e-5 && ges4.objective < 1e-4
@test foldl( ( (x,y) -> x && y ), check_nrm.(ges_nrm_grad))

# MathOptInterface.eval_objective(evaluator, ges1.solution)
# MathOptInterface.eval_objective(evaluator, ges2.solution)
# MathOptInterface.eval_objective_gradient()
