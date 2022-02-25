@testset "PBFGS-PLBFGS" begin 
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

	@test isapprox(norm(ges_mnlp_pbfgs.solution - ones(n), 2), 0, atol=1e-4)
	@test isapprox(norm(ges_adnlp_pbfgs.solution - ones(n), 2), 0, atol=1e-4)
	@test isapprox(norm(ges_mnlp_plbfgs.solution - ones(n), 2), 0, atol=1e-4)
	@test isapprox(norm(ges_adnlp_plbfgs.solution - ones(n), 2), 0, atol=1e-4)
end 