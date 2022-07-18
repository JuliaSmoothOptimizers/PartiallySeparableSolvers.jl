module PartiallySeparableSolvers

using LinearAlgebra, Printf
using ExpressionTreeForge, PartiallySeparableNLPModels, PartitionedStructures
using Krylov, LinearOperators
using NLPModels, ADNLPModels, NLPModelsJuMP
using JuMP, MathOptInterface
using ModelingToolkit
using SolverCore

export PTRUNK # Partitioned update solver

include("tr_cg_part_data.jl")
include("partitioned_methods.jl")

using ..ModPartitionedMethods

end
