# PartiallySeparableSolvers.jl

## Philosophy
PartiallySeparableSolvers.jl implements a partitioned quasi-Newton trust-region solver for unconstrained optimization.

## Compatibility
Julia ≥ 1.6.

## How to install
```
pkg> add https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl
pkg> test PartiallySeparableSolvers
```

## How to use 
See the [tutorial](https://paraynaud.github.io/PartiallySeparableSolvers.jl/dev/tutorial/).

## Dependencies
The module uses [PartiallySeparableNLPModels.jl](https://github.com/JuliaSmoothOptimizers/PartiallySeparableNLPModels.jl) to detect the partially-separable structure and the structure required to partitioned quasi-Newton approximations.