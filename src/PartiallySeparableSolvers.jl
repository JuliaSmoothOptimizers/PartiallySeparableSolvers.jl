module PartiallySeparableSolvers

    using LinearOperators, Krylov, LinearAlgebra

    using NLPModels

    using SolverTools
    using NLPModelsJuMP, JuMP, MathOptInterface

    using Printf
    #=----------------------------------------------------------------------------------------------------------=#
    #Ajout
    using CalculusTreeTools, PartiallySeparableNLPModel
    #=----------------------------------------------------------------------------------------------------------=#
    # ..implementation_expr_tree,  ..PartiallySeparableStructure, ..trait_expr_tree,

    include("quasi_newton.jl")
    include("PartitionnedSolvers.jl")
    export solver_TR_PSR1!, solver_TR_PBFGS!
end
