module ModTrustRegionPartitionedData

using ExpressionTreeForge, PartitionedStructures
# , PartiallySeparableNLPModels
using LinearAlgebra, LinearAlgebra.BLAS, LinearOperators, NLPModels, Krylov
using Printf, SolverCore, SolverTools
using ..Mod_ab_partitioned_data, ..Mod_PQN

export partitionedTrunk

"""
    Counter

Substitute to `NLPModels.Counters` since the partitioned methods don't relie (for now) on `PartiallySeparableNLPModel`s.
It has fields:

* `neval_obj::Int`: count objective evaluations;
* `neval_grad::Int`: count gradient computations;
* `neval_Hprod::Int`: count Hessian-approximation-vector products.
"""
mutable struct Counter
  neval_obj::Int
  neval_grad::Int
  neval_Hprod::Int
end
increase_obj!(c::Counter) = c.neval_obj += 1
increase_grad!(c::Counter) = c.neval_grad += 1
increase_Hv(c::Counter) = c.neval_Hprod += 1

"""
    stats = partitionedTrunk(nlp::AbstractNLPModel, part_data::PartitionedData; max_eval::Int = 10000, max_iter::Int = 10000, start_time::Float64 = time(), max_time::Float64 = 30.0, ϵ::Float64 = 1e-6, name = part_data.name, name_method::String = "Trust-region " * String(name), kwargs...)

Produce a `GenericExecutionStats` for a partitioned quasi-Newton trust-region method.
It requires the partitioned structures of `part_data::PartitionedDate` paired with an `nlp` model.
The counter `nlp.counters` are updated with the informations of `ModTrustRegionPartitionedData.Counter` to ease the definition of a `GenericExecutionStats`.
"""
function partitionedTrunk(
  nlp::AbstractNLPModel,
  part_data::PartitionedData;
  x₀ = get_x(part_data),
  T = eltype(x₀),
  n = get_n(part_data),
  max_eval::Int = 10000,
  max_iter::Int = 10000,
  start_time::Float64 = time(),
  max_time::Float64 = 30.0,
  atol::Real = √eps(eltype(x₀)),
  rtol::Real = √eps(eltype(x₀)),
  ϵ::Float64 = 1e-6,
  η::Float64 = 1e-3,
  η₁::Float64 = 0.75, # > η
  Δ::Float64 = 1.0,
  ϕ::Float64 = 2.0,
  name = part_data.name,
  name_method::String = "Trust-region " * String(name),
  verbose::Int = 0,
  verbose_part_update = false,
  kwargs...,
)
  x = copy(x₀)
  ∇f₀ = evaluate_grad_part_data(part_data, x₀)
  ∇f₀Norm2 = norm(∇f₀, 2)

  cpt = Counter(0, 0, 0)
  verbose > 0 && println("Start: " * name_method)

  iter = 0 # ≈ k
  gₖ = copy(∇f₀)
  ∇fNorm2 = ∇f₀Norm2
  gtmp = similar(gₖ)
  sₖ = similar(x)

  f(x::AbstractVector) = NLPModels.obj(nlp, x)

  fₖ = f(x)

  verbose > 0 && @info log_header(
    [:iter, :f, :dual, :radius, :ratio, :cgstatus],
    [Int, T, T, T, T, Int, String],
    hdr_override = Dict(:f => "f(x)", :dual => "π", :radius => "Δ"),
  )

  cgtol = one(T)  # Must be ≤ 1.
  cgtol = max(rtol, min(T(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))

  ρₖ = -1

  B = LinearOperators.LinearOperator(part_data)

  (verbose > 0) && @info log_row([iter, fₖ, ∇fNorm2, Δ, ρₖ, "initial point"])

  # stop condition
  absolute(n, gₖ, ϵ) = norm(gₖ, 2) > ϵ
  relative(n, gₖ, ϵ, ∇f₀Norm2) = norm(gₖ, 2) > ϵ * ∇fNorm2
  _max_iter(iter, max_iter) = iter < max_iter
  _max_time(start_time) = (time() - start_time) < max_time
  while absolute(n, gₖ, ϵ) &&
          relative(n, gₖ, ϵ, ∇fNorm2) &&
          _max_iter(iter, max_iter) & _max_time(start_time) # stop condition
    iter += 1
    cg_res = Krylov.cg(B, -gₖ, atol = T(atol), rtol = cgtol, radius = T(Δ), itmax = max(2 * n, 50))
    sₖ .= cg_res[1] # the step deduce by cg

    (ρₖ, fₖ₊₁) = compute_ratio(x, fₖ, sₖ, part_data, B, gₖ, f; cpt = cpt) # compute pk

    if ρₖ > η
      x .= x .+ sₖ
      fₖ = fₖ₊₁
      gtmp .= gₖ
      update_nlp!(
        part_data,
        sₖ;
        name = part_data.name,
        verbose = (verbose > 1),
      )
      gₖ .= PartitionedStructures.get_v(get_pg(part_data)) # already build by update_nlp!
      increase_grad!(cpt)
      ∇fNorm2 = norm(gₖ, 2)

    else
      fₖ = fₖ
    end
    # trust region update
    (ρₖ >= η₁ && norm(sₖ, 2) >= 0.8 * Δ) ? Δ = ϕ * Δ : Δ = Δ
    (ρₖ <= η) && (Δ = 1 / ϕ * Δ)

    (verbose > 0) && @info log_row([iter, fₖ, ∇fNorm2, Δ, ρₖ, cg_res[2].status])
  end

  Δt = time() - start_time
  f = f(x)
  g = evaluate_grad_part_data(part_data, x)
  nrm_grad = norm(g, 2)

  if !absolute(n, g, ϵ) || !relative(n, g, ϵ, ∇f₀Norm2)
    status = :first_order
  elseif !_max_iter(iter, max_iter)
    status = :max_eval
  elseif !_max_time(start_time)
    status = :max_time
  else
    status = :unknown
  end

  nlp.counters.neval_obj = cpt.neval_obj
  nlp.counters.neval_grad = cpt.neval_grad
  nlp.counters.neval_hprod = cpt.neval_Hprod

  stats = GenericExecutionStats(
    status,
    nlp,
    solution = x,
    iter = iter,
    dual_feas = nrm_grad,
    objective = f,
    elapsed_time = Δt,
  )
  return stats
end

"""
    ρₖ = compute_ratio(x::AbstractVector{T}, fₖ::T, sₖ::Vector{T}, part_data::PartitionedData, B::AbstractLinearOperator{T}, gₖ::AbstractVector{T}; cpt::Counter = Counter(0, 0, 0))

Compute the ratio between the actual loss and the expected loss using `part_data`, the current point `x` and the step `s`.
`g_k` must be the gradient at `x` and `B` the linear-operator paired to `part_data`.
"""
function compute_ratio(
  x::AbstractVector{T},
  fₖ::T,
  sₖ::Vector{T},
  part_data::PartitionedData,
  B::AbstractLinearOperator{T},
  gₖ::AbstractVector{T},
  f::Function;
  cpt::Counter = Counter(0, 0, 0),
) where {T <: Number}
  mₖ₊₁ = fₖ + dot(gₖ, sₖ) + 1 / 2 * (dot((B * sₖ), sₖ))
  fₖ₊₁ = f(x + sₖ) # the evaluation set partdata.x to x+sₖ
  set_x!(part_data, x) # set x to its real value, mandatoy for the computation of y
  increase_obj!(cpt)
  ρₖ = (fₖ - fₖ₊₁) / (fₖ - mₖ₊₁)
  return (ρₖ, fₖ₊₁)
end

end
