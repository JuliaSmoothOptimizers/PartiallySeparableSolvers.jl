using PartiallySeparableSolvers

using Test, LinearAlgebra
using ExpressionTreeForge, PartiallySeparableNLPModels, PartitionedStructures
using JuMP, MathOptInterface
using ADNLPModels, NLPModels, NLPModelsJuMP
using OptimizationProblems, OptimizationProblems.ADNLPProblems
using ModelingToolkit

using PartiallySeparableSolvers.Mod_partitioned_methods

include("models.jl")
include("pss.jl")
include("partitioned_methods.jl")
include("tr_cg_part_data.jl")
