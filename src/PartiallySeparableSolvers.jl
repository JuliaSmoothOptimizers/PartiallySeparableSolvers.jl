module PartiallySeparableSolvers
    using LinearOperators, Krylov, LinearAlgebra

    using NLPModels, ADNLPModels

    using SolverTools, SolverCore
    using NLPModelsJuMP, JuMP, MathOptInterface

		using ModelingToolkit

    using Printf
    #=----------------------------------------------------------------------------------------------------------=#
    #Ajout
    using CalculusTreeTools, PartiallySeparableNLPModels
    #=----------------------------------------------------------------------------------------------------------=#
    # ..implementation_expr_tree,  ..PartiallySeparableStructure, ..trait_expr_tree,



    include("quasi_newton.jl") #définie les mises à jour BFGS/SR1 ainsi que leurs versions par morceaux
    include("PartitionnedSolvers.jl") #définie les solvers PSR1, PBFGS

    include("partitionned_NLPModel.jl")
    # using ..test_Partitionned_NLPModel

    include("impl_Tr_Cg_Ab.jl") # Définie les solvers LSR1 et LBFGS


    s_a_PBFGS(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PBFGS!(m; kwargs...)
    s_a_PSR1(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PSR1!(m; kwargs...)
    s_a_PBS(m :: T; kwargs... ) where T <: AbstractNLPModel = _solver_TR_PBS!(m; kwargs...)

		PBFGS(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PBFGS!(m; kwargs...)[2]
    PSR1(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PSR1!(m; kwargs...)[2]
    PBS(m :: T; kwargs... ) where T <: AbstractNLPModel = _solver_TR_PBS!(m; kwargs...)[2]
    PTRUNK( nlp :: T; kwargs...) where T <: AbstractNLPModel = my_Part_Trunk(nlp; kwargs...)


    export PBFGS, PSR1, PBS, PTRUNK
    export my_LBFGS, my_LSR1
    export PartionnedNLPModel

    # export test_Partitionned_NLPModel

    # export test_Partitionned_NLPModel.PartionnedNLPModel
    # export test_Partitionned_NLPModel.obj, test_Partitionned_NLPModel.grad!, test_Partitionned_NLPModel.hprod!
end
