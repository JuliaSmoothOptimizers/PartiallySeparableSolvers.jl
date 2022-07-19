# PartiallySeparableSolvers.jl Tutorial


## Solvers exploiting partial separability
PartiallySeparableSolvers.jl defines a solver exploiting automatically the partially-separable structure of $f:\R^n \to \R$
```math
 f(x) = \sum_{i=1}^N f_i (U_i x) , \; f_i : \R^{n_i} \to \R, \; U_i \in \R^{n_i \times n},\; n_i \ll n,
```
where $f$ is the sum of element functions $f_i$.
This solver exploits the partitioned structure of the gradient
```math
\nabla f(x) = \sum_{i=1}^N U_i^\top \nabla f_i (U_i x),
```
and the Hessian 
```math
\nabla^2 f(x) = \sum_{i=1}^N U_i^\top \nabla^2 f_i (U_i x) U_i,
```
to maintain a partitioned a quasi-Newton approximation
```math
B = \sum_{i=1}^N U_i^\top B_{i} U_i,
```
where each $B_{i} \approx \nabla^2 f_i$ has size $n_i \times n_i$.

The usual quasi-Newton updates, such as SR1, BFGS, or their limited-memory variants LSR1 or LBFGS, approximate $\nabla^2 f(x)\approx B$ by low rank updates, which inevitably result in dense approximations.
By relying on element Hessian approximations $B_i$ updated at each iteration, the partitioned updates respect the Hessian sparsity structure and perform updates of rank $\min(\lambda N,n), \; \lambda = 1,2$ depending on the update applied to $B_i$.

See [PartitionedStructures.jl tutorial](https://JuliaSmoothOptimizers.github.io/PartitionedStructures.jl/dev/tutorial/) for more details about partitioned quasi-Newton approximations and how they compare with standard updates.
Some of these partitioned quasi-Newton methods are detailed in the reference below, and some are new (see PartitionedStructures.jl).

#### Reference
* A. Griewank and Ph. L. Toint, [*Partitioned variable metric updates for large structured optimization problems*](https://link.springer.com/article/10.1007/BF01399316), Numerische Mathematik volume, 39, pp. 119--137, 1982.


## Running a partitioned quasi-Newton solver
For now, PartiallySeparableSolvers.jl supports [ADNLPModel](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl)s, which are defined by pure julia code, or [MathOptNLPModel](https://github.com/JuliaSmoothOptimizers/NLPModelsJuMP.jl)s, which are based on [JuMP](https://github.com/jump-dev/JuMP.jl) models.
Regardless of which model is used, the solver returns a [`GenericExecutionStats`](https://juliasmoothoptimizers.github.io/SolverCore.jl/dev/reference/#SolverCore.GenericExecutionStats).

### An `ADNLPModel`
Let's first define the function that we seek to minimize and wrap it into an `ADNLPModel`:
```@example PSSolver
using ADNLPModels

function f(x)
  n = length(x)
  n < 2 && @error("length of x must be >= 2")
  return sum((x[j] + x[j+1])^2 for j=1:n-1)
end 
model = ADNLPModel(f, ones(10))
```

### Running the partitioned quasi-Newton trust-region method
You minimize your model by calling the partitioned-update solver `PTRUNK`:
```@example PSSolver
using PartiallySeparableSolvers
stats = PTRUNK(model)
print(stats)
```

The partitioned quasi-Newton update performed by `PTRUNK` is chosen by the optional argument `name`.
Allowed values include:
- `name=:pbfgs` each $B_i$ is updated by BFGS;
- `name=:psr1` each $B_i$ is updated by SR1;
- `name=:pse` each $B_i$ may be updated by both BFGS and SR1;
- `name=:plbfgs` each $B_i$ is an LBFGS operator;
- `name=:plsr1` each $B_i$ is an LSR1 operator;
- `name=:plse` (by default) each $B_i$ may be an LBFGS or LSR1 operator.

```@example PSSolver
stats = PTRUNK(mathopt_model; name=:pbfgs)
print(stats)
```

See [PartitionedStructures.jl tutorial](https://JuliaSmoothOptimizers.github.io/PartitionedStructures.jl/dev/tutorial/) for more details about partitioned quasi-Newton approximations.