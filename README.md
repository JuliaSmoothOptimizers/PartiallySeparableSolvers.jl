# PartiallySeparableSolvers : Trust-region methods with partitioned quasi-Newton approximations

| **Documentation** | **Linux/macOS/Windows/FreeBSD** | **Coverage** | **DOI** |
|:-----------------:|:-------------------------------:|:------------:|:-------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![build-gh][build-gh-img]][build-gh-url] | [![codecov][codecov-img]][codecov-url] | [![doi][doi-img]][doi-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaSmoothOptimizers.github.io/PartiallySeparableSolvers.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://JuliaSmoothOptimizers.github.io/PartiallySeparableSolvers.jl/dev
[build-gh-img]: https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/workflows/CI/badge.svg?branch=master
[build-gh-url]: https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/actions
[codecov-img]: https://codecov.io/gh/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/branch/master/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl
[doi-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.6984386.svg
[doi-url]: https://zenodo.org/badge/latestdoi/267339899

## How to cite

If you use PartiallySeparableSolvers.jl in your work, please cite using the format given in [CITATION.bib](CITATION.bib).

## Philosophy
PartiallySeparableSolvers.jl implements a partitioned quasi-Newton trust-region solvers for unconstrained optimization.

## Compatibility
Julia ≥ 1.6.

## How to install
```
pkg> add PartiallySeparableSolvers
pkg> test PartiallySeparableSolvers
```

## How to use 
See the [tutorial](https://JuliaSmoothOptimizers.github.io/PartiallySeparableSolvers.jl/dev/tutorial/).

## Dependencies
The module uses [PartiallySeparableNLPModels.jl](https://github.com/JuliaSmoothOptimizers/PartiallySeparableNLPModels.jl) to detect the partially-separable structure and allocate related partitioned structures to perform partitioned quasi-Newton updates.


# Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers), so questions about any of our packages are welcome.
