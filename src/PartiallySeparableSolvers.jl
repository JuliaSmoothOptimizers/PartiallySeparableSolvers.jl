module PartiallySeparableSolvers
    using LinearAlgebra, Printf
		using Krylov, LinearOperators
    using NLPModels, ADNLPModels, NLPModelsJuMP
    using SolverTools, SolverCore
    using JuMP, MathOptInterface
		using ModelingToolkit
    using CalculusTreeTools, PartiallySeparableNLPModels, PartitionedStructures

		include("old_version/_include.jl")
		include("new_version/_include.jl")

		using ..Mod_partitioned_methods

    export PBFGS, PSR1, PBS, PTRUNK
		export PBFGS2, PLBFGS, PUS
    export my_LBFGS, my_LSR1
    export PartionnedNLPModel


    s_a_PBFGS(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PBFGS!(m; kwargs...)
    s_a_PSR1(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PSR1!(m; kwargs...)
    s_a_PBS(m :: T; kwargs... ) where T <: AbstractNLPModel = _solver_TR_PBS!(m; kwargs...)

		PBFGS(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PBFGS!(m; kwargs...)[2]
    PSR1(m :: T;  kwargs... ) where T <: AbstractNLPModel = solver_TR_PSR1!(m; kwargs...)[2]
    PBS(m :: T; kwargs... ) where T <: AbstractNLPModel = _solver_TR_PBS!(m; kwargs...)[2]
    PTRUNK( nlp :: T; kwargs...) where T <: AbstractNLPModel = my_Part_Trunk(nlp; kwargs...)


end
