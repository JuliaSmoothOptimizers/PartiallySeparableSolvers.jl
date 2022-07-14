module PartiallySeparableSolvers
using LinearAlgebra, Printf
using Krylov, LinearOperators
using NLPModels, ADNLPModels, NLPModelsJuMP
using SolverTools, SolverCore
using JuMP, MathOptInterface
using ModelingToolkit
using ExpressionTreeForge, PartiallySeparableNLPModels, PartitionedStructures

include("impl_Tr_Cg_Ab.jl")
include("tr_cg_part_data.jl")
include("partitioned_methods.jl")

using ..Mod_partitioned_methods

export PUS # Partitioned update solver
export my_LBFGS, my_LSR1

end
