@testset "Limited memory" begin
  n = 10
  ros_mnlp = create_Rosenbrock_MathOptModel(n)
  ros_adnlp = Rosenbrock_ADNLPModel(n)

  ges_lbfgs1 = my_LBFGS(ros_mnlp; verbose=false)
  ges_lbfgs2 = my_LBFGS(ros_adnlp; verbose=false)
  ges_lsr11 = my_LSR1(ros_mnlp; verbose=false)
  ges_lsr12 = my_LSR1(ros_adnlp; verbose=false)

  @test ges_lbfgs1.status == :first_order
  @test ges_lbfgs2.status == :first_order
  @test ges_lsr11.status == :first_order
  @test ges_lsr12.status == :first_order
end