using PartiallySeparableSolvers.Mod_ab_partitioned_data
using PartiallySeparableSolvers.Mod_PQN

@testset "PQN structure" begin
  n = 18
  adnlp = ADNLPProblems.dixmaanh(; n)
  obj = ExpressionTreeForge.get_expression_tree(adnlp)

  x = (x -> 2 * x).(ones(n))
  y = rand(n)

  ps_data = build_PartitionedDataTRPQN(obj, n; x0 = x)

  res = show(ps_data)
  @test res == nothing

  objx = evaluate_obj_part_data(ps_data, x)
  @test objx ≈ NLPModels.obj(adnlp, x)

  objy = evaluate_obj_part_data(ps_data, y)
  @test objy ≈ NLPModels.obj(adnlp, y)

  gx = evaluate_grad_part_data(ps_data, x)
  @test NLPModels.grad(adnlp, x) ≈ gx

  gy = evaluate_grad_part_data(ps_data, y)
  @test NLPModels.grad(adnlp, y) ≈ gy

  x = (x -> 2 * x).(ones(n))
  s = (x -> 0.1 * x).(ones(n))

  evaluate_y_part_data!(ps_data, x, s)

  v = ones(n)
  hprod(ps_data, x, v)

  Bk = Matrix(ps_data.pB)

  update_nlp!(ps_data, x, s)
  Bk1 = Matrix(ps_data.pB)
  epv_y = ps_data.py
  PartitionedStructures.build_v!(epv_y)
  _y = PartitionedStructures.get_v(epv_y)
  @test isapprox(norm(Bk1 * s - _y), 0, atol = 1e-10)

  res = product_part_data_x(ps_data, x)
  @test res ≈ Bk1 * x
end

@testset "PartitionedDataTRPQN from ADNLPProblems, update_nlp!(part_data, x, s)" begin
  n = 18
  adnlp = ADNLPProblems.dixmaanh(; n)
  obj = ExpressionTreeForge.get_expression_tree(adnlp)

  x = (x -> 2 * x).(ones(n))
  s = rand(n)

  ps_data_plbfgs = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :plbfgs)
  # ps_data_plbfgs_damped = build_PartitionedDataTRPQN(obj, n;x0=x, name=:plbfgs, damped=true)
  # ps_data_plsr1 = build_PartitionedDataTRPQN(obj, n;x0=x, name=:plsr1)
  ps_data_plse = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :plse)
  ps_data_pbfgs = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pbfgs)
  ps_data_psr1 = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :psr1)
  ps_data_pse = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pse)
  ps_data_pcs = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pcs)

  update_nlp!(ps_data_plbfgs, x, s; verbose = false)
  update_nlp!(ps_data_plse, x, s; verbose = false)
  update_nlp!(ps_data_pbfgs, x, s; verbose = false)
  update_nlp!(ps_data_psr1, x, s; verbose = false)
  update_nlp!(ps_data_pse, x, s; verbose = false)
  update_nlp!(ps_data_pcs, x, s; verbose = false)
  # update_nlp!(ps_data_plbfgs_damped, x, s)
  # update_nlp!(ps_data_plsr1, x, s)

  @test ps_data_plbfgs.py == ps_data_plse.py
  @test ps_data_plbfgs.py == ps_data_pbfgs.py
  @test ps_data_plbfgs.py == ps_data_psr1.py
  @test ps_data_plbfgs.py == ps_data_pse.py
  @test ps_data_plbfgs.py == ps_data_pcs.py
  # @test ps_data_plsr1.py == ps_data_plbfgs.py
  # @test ps_data_plbfgs.py != ps_data_plbfgs_damped.py

  epv_y = ps_data_plbfgs.py
  PartitionedStructures.build_v!(epv_y)
  y = PartitionedStructures.get_v(epv_y)

  partitioned_matrix(nlp) = Matrix(nlp.pB)

  # in the case of the Rosenbrock equation, for the given x,s and induces y, every partitioned update ensure the secant equation.
  @test isapprox(norm(partitioned_matrix(ps_data_plbfgs) * s - y), 0, atol = 1e-10)
  # @test isapprox(norm(partitioned_matrix(ps_data_plbfgs_damped)*s - y_damped), 0, atol=1e-10)
  # @test isapprox(norm(partitioned_matrix(ps_data_plsr1)*s - y), 0, atol=1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_plse) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_pbfgs) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_psr1) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_pse) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_pcs) * s - y), 0, atol = 1e-10)
end

@testset "PartitionedDataTRPQN from PureJuMP, update_nlp!(part_data, s)" begin
  n = 18
  adnlp = PureJuMP.dixmaanh(; n)
  obj = ExpressionTreeForge.get_expression_tree(adnlp)

  x = (x -> 2 * x).(ones(n))
  s = rand(n)

  ps_data_plbfgs = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :plbfgs)
  ps_data_plse = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :plse)
  ps_data_pbfgs = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pbfgs)
  ps_data_psr1 = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :psr1)
  ps_data_pse = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pse)
  ps_data_pcs = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pcs)

  evaluate_grad_part_data!(ps_data_plbfgs)
  evaluate_grad_part_data!(ps_data_plse)
  evaluate_grad_part_data!(ps_data_pbfgs)
  evaluate_grad_part_data!(ps_data_psr1)
  evaluate_grad_part_data!(ps_data_pse)
  evaluate_grad_part_data!(ps_data_pcs)

  update_nlp!(ps_data_plbfgs, s; verbose = false)
  update_nlp!(ps_data_plse, s; verbose = false)
  update_nlp!(ps_data_pbfgs, s; verbose = false)
  update_nlp!(ps_data_psr1, s; verbose = false)
  update_nlp!(ps_data_pse, s; verbose = false)
  update_nlp!(ps_data_pcs, s; verbose = false)

  @test ps_data_plbfgs.py == ps_data_plse.py
  @test ps_data_plbfgs.py == ps_data_pbfgs.py
  @test ps_data_plbfgs.py == ps_data_psr1.py
  @test ps_data_plbfgs.py == ps_data_pse.py
  @test ps_data_plbfgs.py == ps_data_pcs.py

  epv_y = ps_data_plbfgs.py
  PartitionedStructures.build_v!(epv_y)
  y = PartitionedStructures.get_v(epv_y)

  partitioned_matrix(nlp) = Matrix(nlp.pB)

  # in the case of the Rosenbrock equation, for the given x,s and induces y, every partitioned update ensure the secant equation.
  @test isapprox(norm(partitioned_matrix(ps_data_plbfgs) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_plse) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_pbfgs) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_psr1) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_pse) * s - y), 0, atol = 1e-10)
  @test isapprox(norm(partitioned_matrix(ps_data_pcs) * s - y), 0, atol = 1e-10)
end

@testset "methods" begin
  n = 20
  adnlp = ADNLPProblems.arwhead(; n)
  obj = ExpressionTreeForge.get_expression_tree(adnlp)

  x = (x -> 2 * x).(ones(n))
  s = rand(n)

  part_data = build_PartitionedDataTRPQN(obj, n; x0 = x, name = :pbfgs)

  @test get_vec_elt_fun(part_data) == part_data.vec_elt_fun
  @test get_vec_elt_complete_expr_tree(part_data) == part_data.vec_elt_complete_expr_tree
  @test get_element_expr_tree_table(part_data) == part_data.element_expr_tree_table
  @test get_vec_compiled_element_gradients(part_data) == part_data.vec_compiled_element_gradients

  s = rand(n)
  set_s!(part_data, s)
  @test get_s(part_data) == s

  set_n!(part_data, n + 1)
  @test get_n(part_data) == n + 1
  set_n!(part_data, n)

  N = get_N(part_data)
  set_N!(part_data, N + 1)
  @test get_N(part_data) == N + 1
  set_N!(part_data, N)

  v = rand(n)
  set_v!(part_data, v)
  @test Mod_ab_partitioned_data.get_v(part_data) == v

  onesn = ones(n)
  epv = similar(get_pv(part_data))
  epv_from_v!(epv, onesn)

  set_pv!(part_data, onesn)
  @test get_pv(part_data) == epv

  set_pv!(part_data, epv)
  @test get_pv(part_data) == epv

  set_ps!(part_data, onesn)
  @test get_ps(part_data) == epv

  set_pg!(part_data, onesn)
  @test get_pg(part_data) == epv

  set_pg!(part_data, epv)
  @test get_pg(part_data) == epv

  set_py!(part_data, onesn)
  @test get_py(part_data) == epv

  set_phv!(part_data, onesn)
  @test get_phv(part_data) == epv

  epm = epm_from_epv(epv)

  set_pB!(part_data, epm)
  @test get_pB(part_data) == epm
end
