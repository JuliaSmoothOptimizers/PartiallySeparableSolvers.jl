using PartiallySeparableSolvers

using Test, LinearAlgebra, SparseArrays
using CalculusTreeTools, PartiallySeparableNLPModels, PartitionedStructures
using JuMP, MathOptInterface
using ADNLPModels, NLPModels, NLPModelsJuMP
using JSOSolvers, ModelingToolkit

last = true
not_last = true

not_last && include("old_version/_include.jl")
last && include("new_version/_include.jl")
