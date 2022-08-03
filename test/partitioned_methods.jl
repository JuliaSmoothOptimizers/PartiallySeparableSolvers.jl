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
    ex_ = get_expression_tree(ros_mnlp)
    @test ex_ == ex_mathopt    
  end

  @testset "ADNLPModel" begin
    ex_ = get_expression_tree(ros_adnlp)
    @test ex_ == ex_ad    
  end

  @testset "ADNLPModel/MathOptNLPModel" begin
    ex_1 = get_expression_tree(ros_mnlp)
    ex_2 = get_expression_tree(ros_adnlp)

    @test ExpressionTreeForge.evaluate_expr_tree(ex_2, x0) ==
          ExpressionTreeForge.evaluate_expr_tree(ex_1, x0)

    pd_pbfgs1 = build_PartitionedDataTRPQN(ex_1, n; x0)
    pd_pbfgs2 = build_PartitionedDataTRPQN(ex_2, n; x0)
    @test pd_pbfgs1.M == pd_pbfgs2.M
    @test pd_pbfgs1.N == pd_pbfgs2.N
    @test pd_pbfgs1.n == pd_pbfgs2.n
  end
end
