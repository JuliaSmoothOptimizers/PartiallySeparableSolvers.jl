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

    include("quasi_newton.jl") #définie les mises à jour BFGS/SR1 ainsi que leurs versions par morceaux
    include("PartitionnedSolvers.jl") #définie les solvers PSR1, PBFGS
    include("impl_TR_Cg_Ab.jl") # Définie les solvers LSR1 et LBFGS

    export solver_TR_PSR1!, solver_TR_PBFGS!
    export LBFGS, LSR1
end
