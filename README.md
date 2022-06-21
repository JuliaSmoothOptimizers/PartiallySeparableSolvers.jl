# PartiallySeparableSolvers : Trust-region methods with partitioned quasi-Newton approximation

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

## Motivation
The module seeks to minimize the partially separable functions
$$
f(x) = \sum_{=1}^N \hat{f}_i (U_i x), \quad f \in \R^n \to \R, \quad \hat f_i:\R^{n_i} \to \R, \quad U_i \in \R^{n_i \times n}.
$$
$f$ is a sum of element functions $\hat{f}_i$, and usually $n_i \ll n$. $U_i$ is a linear operator, it selects the variables used by $\hat{f}_i$.

PartiallySeparableSolvers.jl define a trust-region method, which compute and approximate the derivatives of $f$ by exploiting the partially separable structure.
The derivatives are partitioned, the gradient 
$$
\nabla f(x) = \sum_{i=1}^N U_i^\top \nabla \hat{f}_i (U_i x),
$$
and the hessian 
$$
\nabla^2 f(x) = \sum_{i=1}^N U_i^\top \nabla^2 \hat{f_i} (U_i x) U_i,
$$
are the sum of the element derivatives $\nabla \hat{f}_i,  \nabla^2\hat{f}_i$.
This structure allows to define a partitioned quasi-Newton approximation of $\nabla^2 f$
$$
B = \sum_{i=1}^N U_i^\top \hat{B}_{i} U_i
$$
such that each $\hat{B}_i \approx \nabla^2 \hat{f}_i$.
Contrary to the BFGS and SR1 updates, a partitioned quasi-Newton approximation keeps the sparsity structure of $\nabla^2 f$, moreover, the rank of update $B$ is proportionnal to $\min(N,n)$ while unstructured quasi-Newton update's rank is 1 or 2.
#### Reference
* A. Griewank and P. Toint, [*Partitioned variable metric updates for large structured optimization problems*](10.1007/BF01399316), Numerische Mathematik volume, 39, pp. 119--137, 1982.


## Content
The module take an [ADNLPModel](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl), which is defined by julia pure code, as a model to be minimized.

Here is a quick example explaining how to use PartiallySeparableSolvers.jl
```julia

using PartiallySeparableSolvers, ADNLPModels

function example(x)
  n = length(x)
	n < 2 && @error("length of x must be >= 2")
	return sum( sum( x[j] for j=1:i)^2 for i=2:n)
end 
start_example(n :: Int) = ones(n)
example_ADNLPModel(n :: Int=100) = ADNLPModel(example, start_example(n), name="Example " * string(n) * " variables")

n = 50 # size of the problem
example_nlp_model = example_ADNLPModel(n) # example model of size n
```

Then, you minimize your problem with `PUS` (partitioned update solver), a trust region method that solves at each iterate a partitioned quasi-Newton subproblem
```julia
PUS(example_nlp_model; name=:plbfgs)
```
The partitioned quasi-Newton operator use in the method is specified by :
- `name=:pbfgs` # each $\hat{B}_i$ is update by BFGS
- `name=:psr1` # each $\hat{B}_i$ is update by SR1
- `name=:pse` # $\hat{B}_i$ may be updated by both BFGS and SR1
- `name=:plbfgs` # each $\hat{B}_i$ is LBFGS operator
- `name=:plsr1` # each $\hat{B}_i$ is LSR1 operator
- `name=:plse` # each $\hat{B}_i$ may be a LBFGS or LSR1 operator

## Dependencies
The module [CalculusTreeTools.jl](https://github.com/paraynaud/CalculusTreeTools.jl) detects automatically the partially separable structure of $f$.
The partitioned quasi-Newton approximation are defined in the module [PartitionedStructures.jl](https://github.com/paraynaud/PartitionedStructures.jl).
All the structures required by the trust-region method are managed by [PartiallySeparableNLPModels.jl](https://github.com/paraynaud/PartiallySeparableNLPModels.jl).


## How to install
```
julia> ]
pkg> add https://github.com/paraynaud/PartitionedStructures.jl, https://github.com/paraynaud/CalculusTreeTools.jl, https://github.com/paraynaud/PartiallySeparableNLPModels.jl, https://github.com/paraynaud/PartiallySeparableSolvers.jl
pkg> test PartiallySeparableSolvers
```