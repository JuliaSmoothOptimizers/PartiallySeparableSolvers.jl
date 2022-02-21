include("quasi_newton.jl") #définie les mises à jour BFGS/SR1 ainsi que leurs versions par morceaux
include("PartitionnedSolvers.jl") #définie les solvers PSR1, PBFGS
include("partitionned_NLPModel.jl")
include("impl_Tr_Cg_Ab.jl") # Définie les solvers LSR1 et LBFGS