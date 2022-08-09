@testset "instantiation new version" begin
  n = 10
  adnlp = ADNLPProblems.arwhead(; n)
  monlp = MathOptNLPModel(PureJuMP.arwhead(; n))

  x0 = monlp.meta.x0

  @testset "ADNLPModel/MathOptNLPModel" begin
    ex_1 = get_expression_tree(monlp)
    ex_2 = get_expression_tree(adnlp)

    @test ExpressionTreeForge.evaluate_expr_tree(ex_2, x0) ==
          ExpressionTreeForge.evaluate_expr_tree(ex_1, x0)

    pd_pbfgs1 = build_PartitionedDataTRPQN(ex_1, n; x0)
    pd_pbfgs2 = build_PartitionedDataTRPQN(ex_2, n; x0)
    @test pd_pbfgs1.M == pd_pbfgs2.M
    @test pd_pbfgs1.N == pd_pbfgs2.N
    @test pd_pbfgs1.n == pd_pbfgs2.n
  end
end
