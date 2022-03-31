@testset "PBFGS-PLBFGS Rosenbrock" begin 
	n = 10
	ros_mnlp = create_Rosenbrock_MathOptModel(n)		
	ros_adnlp = Rosenbrock_ADNLPModel(n)

	ges_mnlp_pbfgs = PBFGS2(ros_mnlp)
	ges_adnlp_pbfgs = PBFGS2(ros_adnlp)

	ges_mnlp_plbfgs = PLBFGS(ros_mnlp)
	ges_adnlp_plbfgs = PLBFGS(ros_adnlp)

	@test ges_mnlp_pbfgs.status == :first_order
	@test ges_adnlp_pbfgs.status == :first_order
	@test ges_mnlp_plbfgs.status == :first_order
	@test ges_adnlp_plbfgs.status == :first_order

	@test isapprox(norm(ges_mnlp_pbfgs.solution - ones(n), 2), 0, atol=1e-2)
	@test isapprox(norm(ges_adnlp_pbfgs.solution - ones(n), 2), 0, atol=1e-2)
	@test isapprox(norm(ges_mnlp_plbfgs.solution - ones(n), 2), 0, atol=1e-2)
	@test isapprox(norm(ges_adnlp_plbfgs.solution - ones(n), 2), 0, atol=1e-2)
end 

@testset "PBFGS-PLBFGS Rosenbrock" begin 
	n = 10
	nlp = MathOptNLPModel(OptimizationProblems.arwhead(n), name="arwhead " * string(n))	

	ges_pbfgs = PBFGS2(nlp)
	ges_plbfgs = PLBFGS(nlp)

	@test ges_pbfgs.status == :first_order
	@test ges_plbfgs.status == :first_order	
end 

@testset "PUS Rosenbrock" begin

	n = 10
	nlp = MathOptNLPModel(OptimizationProblems.arwhead(n), name="arwhead " * string(n))	

	ges_plbfgs = PUS(nlp; name=:plbfgs)
	ges_plsr1 = PUS(nlp; name=:plsr1)
	ges_plse = PUS(nlp; name=:plse)
	ges_psr1 = PUS(nlp; name=:psr1)
	ges_pse = PUS(nlp; name=:pse)
	ges_pbfgs = PUS(nlp; name=:pbfgs)

	@test ges_plbfgs.status == :first_order
	@test ges_plsr1.status == :first_order
	@test ges_plse.status == :first_order
	@test ges_psr1.status == :first_order
	@test ges_pse.status == :first_order
	@test ges_pbfgs.status == :first_order

end