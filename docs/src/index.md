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

# Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers), so questions about any of our packages are welcome.
