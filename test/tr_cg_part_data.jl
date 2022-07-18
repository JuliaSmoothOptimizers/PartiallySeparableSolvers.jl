@testset "PTRUNK Arwhead" begin
  n = 10
  nlp = arwhead(; n)

  ges_plbfgs = PTRUNK(nlp; name = :plbfgs, verbose = false)
  ges_plbfgs_damped = PTRUNK(nlp; name = :plbfgs, damped = true, verbose = false)
  ges_plsr1 = PTRUNK(nlp; name = :plsr1, verbose = false)
  ges_plse = PTRUNK(nlp; name = :plse, verbose = false)
  ges_psr1 = PTRUNK(nlp; name = :psr1, verbose = false)
  ges_pse = PTRUNK(nlp; name = :pse, verbose = false)
  ges_pbfgs = PTRUNK(nlp; name = :pbfgs, verbose = false)

  @test ges_plbfgs.status == :first_order
  @test ges_plbfgs_damped.status == :first_order
  @test ges_plsr1.status == :first_order
  @test ges_plse.status == :first_order
  @test ges_psr1.status == :first_order
  @test ges_pse.status == :first_order
  @test ges_pbfgs.status == :first_order
end
