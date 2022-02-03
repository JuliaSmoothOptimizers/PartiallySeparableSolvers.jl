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

@testset "Hv" begin
	temp_sps = similar(x)
	temp_jump = similar(x)
	H_sps = hess_op!(SPS_nlp, x, temp_sps)
	H_jump = hess_op!(JuMP_mod, x, temp_jump)

	H_sps * v
	@test H_sps * v == H_jump * v
end