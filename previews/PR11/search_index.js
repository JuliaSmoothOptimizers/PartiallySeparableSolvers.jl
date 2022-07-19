var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [PartiallySeparableSolvers, ModPartitionedMethods, ModTrustRegionPartitionedData]","category":"page"},{"location":"reference/#PartiallySeparableSolvers.ModPartitionedMethods.PTRUNK-Tuple{NLPModels.AbstractNLPModel}","page":"Reference","title":"PartiallySeparableSolvers.ModPartitionedMethods.PTRUNK","text":"stats = PTRUNK(nlp::AbstractNLPModel; name = :plse, kwargs...)\n\nPTRUNK (partitioned update solver) return a stats::GenericExecutionStats from a partitioned quasi-Newton trust-region method. Several variants are available via the optional argument name. Each variant approximate the element-Hessians differently, some use dense matrices while others use linear-operators. Variants with dense element-Hessian approximations:\n\nPBFGS with name=:pbfgs, every element-Hessian approximation is updated with BFGS;\nPSR1 with name=:psr1, every element-Hessian approximation is updated with SR1;\nPSE with name=:pse, every element-Hessian approximation is updated with BFGS or SR1 if the curvature condition doesn't hold.\n\nVariants with linear-operator element-Hessian approximations:\n\nPLBFGS with name=:plbfgs, every element-Hessian approximations is a LBFGS operator;\nPLSR1 with name=:plsr1, every element-Hessian approximations is a LSR1 operator;\nPLSE with name=:plse, by default, every element-Hessian approximations is a LBFGS operator as long as the curvature condition holds, otherwise it becomes a LSR1 operator.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PartiallySeparableSolvers.ModPartitionedMethods.get_expr_tree-Union{Tuple{NLPModelsJuMP.MathOptNLPModel}, Tuple{T}} where T<:Number","page":"Reference","title":"PartiallySeparableSolvers.ModPartitionedMethods.get_expr_tree","text":"expr_tree, n, x0 = get_expr_tree(nlp::MathOptNLPModel; x0::Vector{T} = copy(nlp.meta.x0), kwargs...) where {T <: Number}\nexpr_tree, n, x0 = get_expr_tree(adnlp::ADNLPModel; x0::Vector{T} = copy(adnlp.meta.x0), kwargs...) where {T <: Number}\n\nReturn the expr_tree, the size n and the initial point x0 from either a MathOptNLPModel or a ADNLPModel.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PartiallySeparableSolvers.ModTrustRegionPartitionedData.Counter","page":"Reference","title":"PartiallySeparableSolvers.ModTrustRegionPartitionedData.Counter","text":"Counter\n\nSubstitute to NLPModels.Counters since the partitioned methods don't relie (for now) on PartiallySeparableNLPModels. It has fields:\n\nneval_obj::Int: count objective evaluations;\nneval_grad::Int: count gradient computations;\nneval_Hprod::Int: count Hessian-approximation-vector products.\n\n\n\n\n\n","category":"type"},{"location":"reference/#PartiallySeparableSolvers.ModTrustRegionPartitionedData.compute_ratio-Union{Tuple{T}, Tuple{AbstractVector{T}, T, Vector{T}, PartiallySeparableNLPModels.Mod_ab_partitioned_data.PartitionedData, LinearOperators.AbstractLinearOperator{T}, AbstractVector{T}}} where T<:Number","page":"Reference","title":"PartiallySeparableSolvers.ModTrustRegionPartitionedData.compute_ratio","text":"ρₖ = compute_ratio(x::AbstractVector{T}, fₖ::T, sₖ::Vector{T}, part_data::PartiallySeparableNLPModels.PartitionedData, B::AbstractLinearOperator{T}, gₖ::AbstractVector{T}; cpt::Counter = Counter(0, 0, 0))\n\nCompute the ratio between the actual loss and the expected loss using part_data, the current point x and the step s. g_k must be the gradient at x and B the linear-operator paired to part_data.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PartiallySeparableSolvers.ModTrustRegionPartitionedData.partitionedTrunk-Tuple{NLPModels.AbstractNLPModel, PartiallySeparableNLPModels.Mod_ab_partitioned_data.PartitionedData}","page":"Reference","title":"PartiallySeparableSolvers.ModTrustRegionPartitionedData.partitionedTrunk","text":"stats = partitionedTrunk(nlp::AbstractNLPModel, part_data::PartiallySeparableNLPModels.PartitionedData; max_eval::Int = 10000, max_iter::Int = 10000, start_time::Float64 = time(), max_time::Float64 = 30.0, ϵ::Float64 = 1e-6, name = part_data.name, name_method::String = \"Trust-region \" * String(name), kwargs...)\n\nProduce a GenericExecutionStats for a partitioned quasi-Newton trust-region method. It requires the partitioned structures of part_data::PartitionedDate paired with an nlp model. The counter nlp.counters are updated with the informations of ModTrustRegionPartitionedData.Counter to ease the definition of a GenericExecutionStats.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PartiallySeparableSolvers.ModTrustRegionPartitionedData.partitionedTrunkCore-Tuple{PartiallySeparableNLPModels.Mod_ab_partitioned_data.PartitionedData}","page":"Reference","title":"PartiallySeparableSolvers.ModTrustRegionPartitionedData.partitionedTrunkCore","text":"(x, iter) = partitionedTrunkCore(part_data::PartiallySeparableNLPModels.PartitionedData; x::AbstractVector = copy(get_x(part_data)), n::Int = get_n(part_data), max_eval::Int = 10000, max_iter::Int = 10000, max_time::Float64 = 30.0, atol::Real = √eps(eltype(x)), rtol::Real = √eps(eltype(x)), start_time::Float64 = time(), η::Float64 = 1e-3, η₁::Float64 = 0.75, # > η Δ::Float64 = 1.0, ϵ::Float64 = 1e-6, ϕ::Float64 = 2.0, ∇f₀::AbstractVector = PartiallySeparableNLPModels.evaluate_grad_part_data(part_data, x), cpt::Counter = Counter(0, 0, 0), iter_print::Int64 = Int(floor(max_iter / 100)), T = eltype(x), verbose = true, kwargs...)\n\nPartitioned quasi-Newton trust-region method apply on the partitioned structures of part_data. The method return the point x and the number of iterations performed before it reaches the stopping criterias.\n\n\n\n\n\n","category":"method"},{"location":"#PartiallySeparableSolvers.jl","page":"Home","title":"PartiallySeparableSolvers.jl","text":"","category":"section"},{"location":"#Philosophy","page":"Home","title":"Philosophy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PartiallySeparableSolvers.jl implements a partitioned quasi-Newton trust-region solver for unconstrained optimization.","category":"page"},{"location":"#Compatibility","page":"Home","title":"Compatibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Julia ≥ 1.6.","category":"page"},{"location":"#How-to-install","page":"Home","title":"How to install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/JuliaSmoothOptimizers/PartiallySeparableSolvers.jl\npkg> test PartiallySeparableSolvers","category":"page"},{"location":"#How-to-use","page":"Home","title":"How to use","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See the tutorial.","category":"page"},{"location":"#Dependencies","page":"Home","title":"Dependencies","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The module uses PartiallySeparableNLPModels.jl to detect the partially-separable structure and the structure required to partitioned quasi-Newton approximations.","category":"page"},{"location":"tutorial/#PartiallySeparableSolvers.jl-Tutorial","page":"Tutorial","title":"PartiallySeparableSolvers.jl Tutorial","text":"","category":"section"},{"location":"tutorial/#Solvers-exploiting-partial-separability","page":"Tutorial","title":"Solvers exploiting partial separability","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"PartiallySeparableSolvers.jl defines a solver exploiting automatically the partially-separable structure of fR^n to R","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":" f(x) = sum_i=1^N f_i (U_i x)   f_i  R^n_i to R  U_i in R^n_i times n n_i ll n","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where f is the sum of element functions f_i. This solver exploits the partitioned structure of the gradient","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nabla f(x) = sum_i=1^N U_i^top nabla f_i (U_i x)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"and the Hessian ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nabla^2 f(x) = sum_i=1^N U_i^top nabla^2 f_i (U_i x) U_i","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"to maintain a partitioned a quasi-Newton approximation","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"B = sum_i=1^N U_i^top B_i U_i","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where each B_i approx nabla^2 f_i has size n_i times n_i.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The usual quasi-Newton updates, such as SR1, BFGS, or their limited-memory variants LSR1 or LBFGS approximate nabla^2 f(x)approx B by low rank updates producing inevitably dense approximations. By relying on element Hessian approximations B_i, updated at each iteration, the partitioned updates respect the Hessian sparsity structure and perform updates of rank min(lambda Nn)  lambda = 12 depending on the update applied to B_i.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"See PartitionedStructures.jl tutorial to get more details about partitioned quasi-Newton approximations and how they compare with standard updates. Some of these partitioned quasi-Newton methods are detailed in the reference below, and some are new (see PartitionedStructures.jl).","category":"page"},{"location":"tutorial/#Reference","page":"Tutorial","title":"Reference","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"A. Griewank and P. Toint, Partitioned variable metric updates for large structured optimization problems, Numerische Mathematik volume, 39, pp. 119–137, 1982.","category":"page"},{"location":"tutorial/#Running-a-partitioned-quasi-Newton-solver","page":"Tutorial","title":"Running a partitioned quasi-Newton solver","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For now, PartiallySeparableSolvers.jl supports ADNLPModels, which are defined by pure julia code, or MathOptNLPModels based on JuMP models. Regardless of which model is used, the solver returns a GenericExecutionStats.","category":"page"},{"location":"tutorial/#An-ADNLPModel","page":"Tutorial","title":"An ADNLPModel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Let's first define the function that we seek to minimize and wrap it into an ADNLPModel:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using ADNLPModels\n\nfunction example(x)\n  n = length(x)\n  n < 2 && @error(\"length of x must be >= 2\")\n  return sum((x[j] + x[j+1])^2 for j=1:n-1)\nend \nstart_example(n :: Int) = ones(n)\nexample_ADNLPModel(n :: Int=100) = ADNLPModel(example, start_example(n), name=\"Example \" * string(n) * \" variables\")\n\nn = 10 # size of the problem\nadnlp_model = example_ADNLPModel(n)","category":"page"},{"location":"tutorial/#Running-the-partitioned-quasi-Newton-trust-region-method","page":"Tutorial","title":"Running the partitioned quasi-Newton trust-region method","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"You minimize your model by calling the partitioned-update solver PTRUNK:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using PartiallySeparableSolvers\nstats = PTRUNK(adnlp_model)\nprint(stats)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The partitioned quasi-Newton update performed by PTRUNK is chosen by the optional argument name. Allowed values include:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"name=:pbfgs each B_i is updated by BFGS;\nname=:psr1 each B_i is updated by SR1;\nname=:pse each B_i may be updated by both BFGS and SR1;\nname=:plbfgs each B_i is an LBFGS operator;\nname=:plsr1 each B_i is an LSR1 operator;\nname=:plse (by default) each B_i may be an LBFGS or LSR1 operator.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"stats = PTRUNK(mathopt_model; name=:pbfgs)\nprint(stats)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"See PartitionedStructures.jl tutorial for more details about partitioned quasi-Newton approximations.","category":"page"}]
}
