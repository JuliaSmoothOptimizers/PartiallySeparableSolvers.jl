# PartiallySeparableSolvers : Trust-region methods with partitioned quasi-Newton approximations

| **Documentation** | **Linux/macOS/Windows/FreeBSD** | **Coverage** | **DOI** |
|:-----------------:|:-------------------------------:|:------------:|:-------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![build-gh][build-gh-img]][build-gh-url] [![build-cirrus][build-cirrus-img]][build-cirrus-url] | [![codecov][codecov-img]][codecov-url] | [![doi][doi-img]][doi-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://paraynaud.github.io/PartiallySeparableSolvers.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://paraynaud.github.io/PartiallySeparableSolvers.jl/dev
[build-gh-img]: https://github.com/paraynaud/PartiallySeparableSolvers.jl/workflows/CI/badge.svg?branch=master
[build-gh-url]: https://github.com/paraynaud/PartiallySeparableSolvers.jl/actions
[build-cirrus-img]: https://img.shields.io/cirrus/github/paraynaud/PartiallySeparableSolvers.jl?logo=Cirrus%20CI
[build-cirrus-url]: https://cirrus-ci.com/github/paraynaud/PartiallySeparableSolvers.jl
[codecov-img]: https://codecov.io/gh/paraynaud/PartiallySeparableSolvers.jl/branch/master/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/paraynaud/PartiallySeparableSolvers.jl
[doi-img]: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.822073-blue.svg
[doi-url]: https://doi.org/10.5281/zenodo.822073


## Philosophy
PartiallySeparableSolvers.jl implements a partitioned quasi-Newton trust-region solvers for unconstrained optimization.

## Compatibility
Julia ≥ 1.6.

## How to install
```
pkg> add https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl
pkg> test PartiallySeparableSolvers
```

## How to use 
See the [tutorial](https://JuliaSmoothOptimizers.github.io/PartiallySeparableSolvers.jl/dev/tutorial/).

## Dependencies
The module uses [PartiallySeparableNLPModels.jl](https://github.com/JuliaSmoothOptimizers/PartiallySeparableNLPModels.jl) to detect the partially-separable structure and the structure required for partitioned quasi-Newton approximations.