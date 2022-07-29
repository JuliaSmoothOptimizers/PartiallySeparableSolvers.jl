# PartiallySeparableSolvers : Trust-region methods with partitioned quasi-Newton approximations

| **Documentation** | **Linux/macOS/Windows/FreeBSD** | **Coverage** | **DOI** |
|:-----------------:|:-------------------------------:|:------------:|:-------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![build-gh][build-gh-img]][build-gh-url] [![build-cirrus][build-cirrus-img]][build-cirrus-url] | [![codecov][codecov-img]][codecov-url] | [![doi][doi-img]][doi-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaSmoothOptimizers.github.io/PartiallySeparableSolvers.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://JuliaSmoothOptimizers.github.io/PartiallySeparableSolvers.jl/dev
[build-gh-img]: https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/workflows/CI/badge.svg?branch=master
[build-gh-url]: https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/actions
[build-cirrus-img]: https://img.shields.io/cirrus/github/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl?logo=Cirrus%20CI
[build-cirrus-url]: https://cirrus-ci.com/github/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl
[codecov-img]: https://codecov.io/gh/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/branch/master/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl
[doi-img]: https://zenodo.org/badge/267339899.svg
[doi-url]: https://zenodo.org/badge/latestdoi/267339899

## How to cite

If you use PartiallySeparableSolvers.jl in your work, please cite using the format given in [CITATION.bib](CITATION.bib).

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
