module PartiallySeparableSolvers

using LinearAlgebra, Printf
using ExpressionTreeForge, PartiallySeparableNLPModels, PartitionedStructures
using Krylov, LinearOperators
using NLPModels, ADNLPModels, NLPModelsJuMP
using SolverCore

export PTRUNK # Partitioned update solver

include("M_partitioned_data.jl")
include("partitioned_data_pqn.jl")

include("tr_cg_part_data.jl")
include("partitioned_methods.jl")

using ..ModPartitionedMethods

end
