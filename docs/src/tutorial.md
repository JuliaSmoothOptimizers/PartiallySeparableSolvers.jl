# PartiallySeparableSolvers.jl Tutorial


## Solvers exploiting partial separability
PartiallySeparableSolvers.jl defines a solver exploiting automatically the partially-separable structure of $f:\R^n \to \R$
```math
 f(x) = \sum_{i=1}^N f_i (U_i x) , \; f_i : \R^{n_i} \to \R, \; U_i \in \R^{n_i \times n},\; n_i \ll n,
```
where $f$ is the sum of element functions $f_i$.
This solver exploit the partitioned structure of the gradient
```math
\nabla f(x) = \sum_{i=1}^N U_i^\top \nabla f_i (U_i x),
```
and the Hessian 
```math
\nabla^2 f(x) = \sum_{i=1}^N U_i^\top \nabla^2 f_i (U_i x) U_i,
```
to produce partitioned a quasi-Newton approximation
```math
B = \sum_{i=1}^N U_i^\top B_{i} U_i,
```
where each $B_{i} \approx \nabla^2 f_i$.

The usual quasi-Newton updates SR1, BFGS or limited-memory variant LSR1, LBFGS approximate $\nabla^2 f(x)$ by low rank update producing inevitably dense approximations.
By relying on element Hessian approximations $B_i$, updated at each iterate, the partitioned updates keep the Hessian sparse structure and perform updates of rank $\min(\lambda N,n), \; \lambda = 1,2$ depending on the update applied to $B_i$.

See [PartitionedStructures.jl tutorial](https://JuliaSmoothOptimizers.github.io/PartitionedStructures.jl/dev/tutorial/) to get more details on how partitioned quasi-Newton approximate $\nabla^2 f$ compare to BFGS.

#### Reference
* A. Griewank and P. Toint, [*Partitioned variable metric updates for large structured optimization problems*](10.1007/BF01399316), Numerische Mathematik volume, 39, pp. 119--137, 1982.


## Run a partitioned quasi-Newton solver
For now, PartiallySeparableSolvers.jl supports either an [ADNLPModel](https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl), which is defined by julia pure code, or a [MathOptNLPModel](https://github.com/JuliaSmoothOptimizers/NLPModelsJuMP.jl) based on a [JuMP](https://github.com/jump-dev/JuMP.jl) model.
Wether the model used, a solver run will return a [`GenericExecutionStats`](https://juliasmoothoptimizers.github.io/SolverCore.jl/dev/reference/#SolverCore.GenericExecutionStats).

### An `ADNLPModel`
You first define the function that you see to minimize and wrap it an `ADNLPModel`:
```@example PSSolver
using ADNLPModels

function example(x)
  n = length(x)
  n < 2 && @error("length of x must be >= 2")
  return sum((x[j] + x[j+1])^2 for j=1:n-1)
end 
start_example(n :: Int) = ones(n)
example_ADNLPModel(n :: Int=100) = ADNLPModel(example, start_example(n), name="Example " * string(n) * " variables")

n = 10 # size of the problem
adnlp_model = example_ADNLPModel(n)
```
### An `MathOptNLPModel`
The same way, you can define a `MathOptNLPModel`:
```@example PSSolver
using JuMP, NLPModelsJuMP

model_jump = Model()
@variable(model_jump, x[1:n])
@NLobjective(model_jump, Min, sum((x[j] + x[j+1])^2 for j = 1:n-1))
variables = JuMP.all_variables(model_jump)
JuMP.set_start_value.(variables, ones(n))
mathopt_model = MathOptNLPModel(model_jump, name = "Example " * string(n))
```

### Run a partitioned quasi-Newton trust-region method
You minimize your model by calling the partitioned update solver `PUS`:
```@example PSSolver
using PartiallySeparableSolvers
PUS(adnlp_model)
```

```@example PSSolver
PUS(mathopt_model)
```

The partitioned quasi-Newton update performed by `PUS` is chosen by the optional argument `name`.
By choosing:
- `name=:pbfgs` each $B_i$ is updated by BFGS;
- `name=:psr1` each $B_i$ is updated by SR1;
- `name=:pse` each $B_i$ may be updated by both BFGS and SR1;
- `name=:plbfgs` each $B_i$ is LBFGS operator;
- `name=:plsr1` each $B_i$ is LSR1 operator;
- `name=:plse` (by default) each $B_i$ may be a LBFGS or LSR1 operator.

```@example PSSolver
PUS(mathopt_model; name=:pbfgs) # by default PLSE
```

See [PartitionedStructures.jl tutorial](https://JuliaSmoothOptimizers.github.io/PartitionedStructures.jl/dev/tutorial/) to get more details about partitioned quasi-Newton approximations.