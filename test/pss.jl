@testset "element function decomposition" begin
  n = 10
  
  nlp = arwhead(;n)
  
  (ex, n, x0) = get_expr_tree(nlp)
  part_data_pbfgs = build_PartitionedDataTRPQN(ex, n; x0 = x0)

  pbfgs_index_element_tree = part_data_pbfgs.index_element_tree
  index_element_tree =
    [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
  @test pbfgs_index_element_tree == index_element_tree
end
