using PartiallySeparableSolvers

using Test, LinearAlgebra, SparseArrays
using ExpressionTreeForge, PartiallySeparableNLPModels, PartitionedStructures
using JuMP, MathOptInterface
using ADNLPModels, NLPModels, NLPModelsJuMP
using JSOSolvers, ModelingToolkit


include("new_version/_include.jl")
