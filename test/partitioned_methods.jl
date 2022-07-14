using PartiallySeparableSolvers.Mod_partitioned_methods

@testset "instantiation new version" begin
  n = 10
  ros_mnlp = create_Rosenbrock_MathOptModel(n)
  ros_adnlp = Rosenbrock_ADNLPModel(n)

  model = ros_mnlp.eval.m
  evaluator = JuMP.NLPEvaluator(model)
  MathOptInterface.initialize(evaluator, [:ExprGraph])
  obj_Expr = MathOptInterface.objective_expr(evaluator)::Expr
  ex_mathopt = ExpressionTreeForge.transform_to_expr_tree(obj_Expr)

  ModelingToolkit.@variables x[1:n]
  fun = ros_adnlp.f(x)
  ex_ad = ExpressionTreeForge.transform_to_expr_tree(fun)
  x0 = create_initial_point_Rosenbrock(n)

  @testset "MathOptNLPModel" begin
    ex_, n_, x0_ = get_expr_tree(ros_mnlp)
    @test ex_ == ex_mathopt
    @test n_ == n
    @test x0_ == x0
  end

  # n = ros_adnlp.meta.nvar
  # ModelingToolkit.@variables x[1:n]
  # fun = ros_adnlp.f(x)
  # ex = ExpressionTreeForge.transform_to_expr_tree(fun)
  @testset "ADNLPModel" begin
    ex_, n_, x0_ = get_expr_tree(ros_adnlp)
    @test ex_ == ex_ad
    @test n_ == n
    @test x0_ == x0
  end

  @testset "ADNLPModel/MathOptNLPModel" begin
    ex_1, n_1, x0_1 = get_expr_tree(ros_mnlp)
    ex_2, n_2, x0_2 = get_expr_tree(ros_adnlp)

    @test ExpressionTreeForge.evaluate_expr_tree(ex_2, x0) ==
          ExpressionTreeForge.evaluate_expr_tree(ex_1, x0)
    @test n_1 == n_2
    @test x0_1 == x0_2

    pd_pbfgs1 = build_PartitionedDataTRPQN(ex_1, n; x0)
    pd_pbfgs2 = build_PartitionedDataTRPQN(ex_2, n; x0)
    @test pd_pbfgs1.M == pd_pbfgs2.M
    @test pd_pbfgs1.N == pd_pbfgs2.N
    @test pd_pbfgs1.n == pd_pbfgs2.n
  end
end
