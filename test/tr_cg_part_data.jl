@testset "PTRUNK Arwhead" begin
  n = 10
  nlp = ADNLPProblems.arwhead(; n)

  ges_plbfgs = PTRUNK(nlp; name = :plbfgs)
  ges_plbfgs_damped = PTRUNK(nlp; name = :plbfgs, damped = true)
  ges_plsr1 = PTRUNK(nlp; name = :plsr1)
  ges_plse = PTRUNK(nlp; name = :plse)
  ges_psr1 = PTRUNK(nlp; name = :psr1)
  ges_pse = PTRUNK(nlp; name = :pse)
  ges_pbfgs = PTRUNK(nlp; name = :pbfgs)

  @test ges_plbfgs.status == :first_order
  @test ges_plbfgs_damped.status == :first_order
  @test ges_plsr1.status == :first_order
  @test ges_plse.status == :first_order
  @test ges_psr1.status == :first_order
  @test ges_pse.status == :first_order
  @test ges_pbfgs.status == :first_order
end
