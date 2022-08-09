using PartiallySeparableSolvers

using Test, LinearAlgebra
using ExpressionTreeForge, PartiallySeparableNLPModels
using ADNLPModels, NLPModels, NLPModelsJuMP
using OptimizationProblems, OptimizationProblems.ADNLPProblems, OptimizationProblems.PureJuMP

using PartiallySeparableSolvers.ModPartitionedMethods

include("pss.jl")
include("partitioned_methods.jl")
include("tr_cg_part_data.jl")
