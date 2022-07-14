using LinearOperators, LinearAlgebra, Krylov, NLPModels
using Printf
using SolverTools

function compute_ratio(
  x::AbstractVector{Y},
  f_x::Y,
  s::Vector{Y},
  nlp::AbstractNLPModel,
  B::AbstractLinearOperator{Y},
  g::AbstractVector{Y},
) where {Y <: Number}
  quad_model_s = f_x + dot(g, s) + 1 / 2 * dot((B * s), s)::Y
  f_next_x = NLPModels.obj(nlp, x + s)::Y
  actual_decrease = f_x - f_next_x::Y
  expected_decrease = f_x - quad_model_s::Y
  return (actual_decrease / expected_decrease, f_next_x)::Tuple{Y, Y}
end

"""
    upgrade_TR_LO(ρₖ, xₖ , sₖ, gₖ, Bₖ, nlp, Δ)

Update the linear operator `Bₖ` and the trust-region radius Δ ratio depending `ρₖ` computed beforhand from the step `sₖ`.
"""
function upgrade_TR_LO!(
  pk::Float64, # value of the ratio
  x_k::AbstractVector{T}, # actual point
  s_k::AbstractVector{T}, # point found at the iteration k
  g_k::AbstractVector{T}, 
  y_k::AbstractVector{T}, 
  B_k::AbstractLinearOperator{T}, # quasi-Newton opertor
  nlp::AbstractNLPModel,
  Δ::Float64; # radius
  η::Float64 = 1e-3,
  η1::Float64 = 0.75
) where {T <: Number}  
  if pk > η
    x_k .= x_k .+ s_k
    y_k .= .- g_k
    NLPModels.grad!(nlp, x_k, g_k)
    y_k .+= g_k
    push!(B_k, s_k, y_k)
  else
  end
  if pk >= η1 # update ∆
    if norm(s_k, 2) < 0.8 * Δ
      Δ = Δ
    else
      Δ = 2 * Δ
    end
  elseif pk <= η
    Δ = 1 / 2 * Δ
  end
  return Δ
end

"""
    solver_TR_CG_Ab_NLP_LO(nlp, B, x0)

A quasi-Newton trust-region solver with conjugate gradient where `B` a limited-memory Hessian approximation.
The method return the final point and the number of iteration performed before it reaches the stopping condition.
"""
function solver_TR_CG_Ab_NLP_LO(
  nlp::AbstractNLPModel,
  B::AbstractLinearOperator{T};
  x::AbstractVector = copy(nlp.meta.x0),
  max_eval::Int = 10000,
  atol::Real = √eps(eltype(x)),
  rtol::Real = √eps(eltype(x)),
  verbose=true,
  kwargs...,
) where {T <: Number}
  T2 = eltype(x)
  (η, η1, Δ, ϵ, iter) = (1e-3, 0.75, 1.0, 10^-6, 0)
  n = nlp.meta.nvar

  g = Vector{T}(undef, length(x))
  ∇f₀ = Vector{T}(undef, length(x))
  sk = similar(g)
  yk = similar(g)
  NLPModels.grad!(nlp, x, ∇f₀)
  g .= ∇f₀
  ∇fNorm2 = nrm2(n, ∇f₀)

  f_xk = NLPModels.obj(nlp, x)
  verbose && (@printf "%3d %8.1e %7.1e %7.1e  \n" iter f_xk norm(g, 2) Δ)

  cgtol = one(T2)  # Must be ≤ 1.
  cgtol = max(rtol, min(T2(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))

  while (nrm2(n, g) > ϵ) && (nrm2(n, g) > ϵ * ∇fNorm2) && iter < max_eval  # stop condition
    iter = iter + 1

    cg_res = Krylov.cg(B, -g, atol = T2(atol), rtol = cgtol, radius = Δ, itmax = max(2 * n, 50))
    sk .= cg_res[1]  # result of the linear system solved by Krylov.cg

    (pk, f_temp) = compute_ratio(x, f_xk, sk, nlp, B, g; η, η1) # we compute the ratio

    Δ = upgrade_TR_LO!(pk, x, sk, g, yk, B, nlp, Δ) # upgrade x,g,B,∆

    if pk > η
      f_xk = f_temp
    end
    verbose && (mod(iter, 500) == 0) && (@printf "%3d %8.1e %7.1e %7.1e  \n" iter f_xk norm(g, 2) Δ)
    
  end

  verbose && (@printf "%3d %8.1e %7.1e %7.1e  \n" iter f_xk norm(g, 2) Δ)

  return (x, iter)
end

#=
MISE EN FORME DES RESULTATS
=#

function solver_TR_CG_Ab_NLP_LO_res(
  nlp::AbstractNLPModel,
  B::AbstractLinearOperator{T};
  max_eval::Int = 10000,
  kwargs...,
) where {T <: Number}
  Δt = @timed ((x_final, iter) = solver_TR_CG_Ab_NLP_LO(nlp, B; max_eval = max_eval, kwargs...))
  x_init = nlp.meta.x0
  nrm_grad = norm(NLPModels.grad(nlp, x_final), 2)
  nrm_grad_init = norm(NLPModels.grad(nlp, x_init), 2)
  if nrm_grad < nrm_grad_init * 1e-6 || nrm_grad < 1e-6
    status = :first_order
  elseif iter >= max_eval
    status = :max_eval
  else
    status = :unknown
  end
  return GenericExecutionStats(
    status,
    nlp,
    solution = x_final,
    iter = iter,
    dual_feas = nrm_grad,
    objective = NLPModels.obj(nlp, x_final),
    elapsed_time = Δt[2],
  )
end

function my_LSR1(
  nlp::AbstractNLPModel;
  x::AbstractVector = copy(nlp.meta.x0),
  kwargs...,
)
  T = eltype(x)
  B = LSR1Operator(nlp.meta.nvar, scaling = true)::LSR1Operator{T} #scaling=true
  println("\n\t LSR1 trust-region method")
  return solver_TR_CG_Ab_NLP_LO_res(nlp, B; x = x, kwargs...)
end

#=
FIN DE LA PARTIE L-SR1
-------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------
Création du LinearOperators L-BFGS, et mise en forme des résultats
Une première partie sans paramètres nommés
=#

function my_LBFGS(
  nlp::AbstractNLPModel;
  x::AbstractVector = copy(nlp.meta.x0),
  kwargs...,
)
  T = eltype(x)
  B = LBFGSOperator(nlp.meta.nvar, scaling = true)::LBFGSOperator{T} #scaling=true
  println("\n\t LBFGS trust-region method")
  return solver_TR_CG_Ab_NLP_LO_res(nlp, B; x = x, kwargs...)
end