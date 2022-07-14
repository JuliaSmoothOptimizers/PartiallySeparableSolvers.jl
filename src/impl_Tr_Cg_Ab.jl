using LinearOperators, LinearAlgebra, Krylov, NLPModels
using Printf

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
  actual_loss = f_x - f_next_x::Y
  expected_loss = f_x - quad_model_s::Y
  return (actual_loss / expected_loss, f_next_x)::Tuple{Y, Y}
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
  max_iter::Int = 10000,
  start_time::Float64 = time(),
  max_time::Float64 = 30.0,
  atol::Real = √eps(eltype(x)),
  rtol::Real = √eps(eltype(x)),
  verbose=true,
  kwargs...,
) where {T <: Number}
  (η, η1, Δ, ϵ, iter) = (1e-3, 0.75, 1.0, 10^-6, 0)
  n = nlp.meta.nvar

  g = Vector{T}(undef, length(x))
  ∇f₀ = Vector{T}(undef, length(x))
  sk = similar(g)
  yk = similar(g)
  NLPModels.grad!(nlp, x, ∇f₀)
  g .= ∇f₀
  ∇f₀Norm2 = norm(∇f₀, 2)

  f_xk = NLPModels.obj(nlp, x)
  verbose && (@printf "%3d %8.1e %7.1e %7.1e  \n" iter f_xk ∇f₀Norm2 Δ)

  cgtol = one(T)  # Must be ≤ 1.
  cgtol = max(rtol, min(T(0.1), 9 * cgtol / 10, sqrt(∇f₀Norm2)))

  absolute(n, gₖ, ϵ) = norm(gₖ, 2) > ϵ
  relative(n, gₖ, ϵ, ∇f₀Norm2) = norm(gₖ, 2) > ϵ * ∇f₀Norm2
  _max_iter(iter, max_iter) = iter < max_iter
  _max_time(start_time) = (time() - start_time) < max_time
  while absolute(n, g, ϵ) &&
          relative(n, g, ϵ, ∇f₀Norm2) &&
          _max_iter(iter, max_iter) & _max_time(start_time)
          
    iter = iter + 1

    cg_res = Krylov.cg(B, -g, atol = T(atol), rtol = cgtol, radius = Δ, itmax = max(2 * n, 50))
    sk .= cg_res[1] # the step deduce by cg

    (pk, f_temp) = compute_ratio(x, f_xk, sk, nlp, B, g) # compute pk

    Δ = upgrade_TR_LO!(pk, x, sk, g, yk, B, nlp, Δ; η, η1) # upgrade x, g, B, ∆
    (pk > η) && (f_xk = f_temp)
        
    verbose && (mod(iter, 50) == 0) && (@printf "%3d %8.1e %7.1e %7.1e  \n" iter f_xk norm(g, 2) Δ)
    
  end

  verbose && (@printf "%3d %8.1e %7.1e %7.1e  \n" iter f_xk norm(g, 2) Δ)

  return (x, iter)
end

"""
    ges = solver_TR_CG_Ab_NLP_LO_ges(nlp::AbstractNLPModel, B::AbstractLinearOperator{T}; max_eval::Int = 10000, max_iter::Int = 10000, start_time::Float64 = time(), max_time::Float64 = 30.0, ϵ::Float64 = 1e-6, kwargs..., ) where {T <: Number}

Return a `GenericExecutionStats` from the quasi-Newton trust-region method `solver_TR_CG_Ab_NLP_LO`.
"""
function solver_TR_CG_Ab_NLP_LO_ges(
  nlp::AbstractNLPModel,
  B::AbstractLinearOperator{T};
  max_eval::Int = 10000,
  max_iter::Int = 10000,
  start_time::Float64 = time(),
  max_time::Float64 = 30.0,
  ϵ::Float64 = 1e-6,
  kwargs...,
) where {T <: Number}
  x_init = nlp.meta.x0
  n = length(x_init)
  ∇f₀ = NLPModels.grad(nlp, x_init)
  ∇f₀Norm2 = norm(∇f₀, 2)

  (x_final, iter) = solver_TR_CG_Ab_NLP_LO(nlp, B; max_eval = max_eval, kwargs...)
  Δt = time() - start_time
  
  f = NLPModels.obj(nlp, x_final)
  g = NLPModels.grad(nlp, x_final)
  ∇fNorm2 = norm(g, 2)

  absolute(n, ∇fNorm2, ϵ) = ∇fNorm2 < ϵ
  relative(n, ∇fNorm2, ϵ, ∇f₀Norm2) = ∇fNorm2 < ϵ * ∇f₀Norm2
  _max_iter(iter, max_iter) = iter >= max_iter
  _max_time(start_time) = (time() - start_time) >= max_time

  if absolute(n, ∇fNorm2, ϵ) || relative(n, ∇fNorm2, ϵ, ∇f₀Norm2)
    status = :first_order
    println("stationnary point ✅")
  elseif _max_iter(iter, max_iter)
    status = :max_eval
    println("Max eval ❌")
  elseif _max_time(start_time)
    status = :max_time
    println("Max time ❌")
  else
    status = :unknown
    println("Unknown ❌")
  end
  ges = GenericExecutionStats(
    status,
    nlp,
    solution = x_final,
    iter = iter,
    dual_feas = ∇fNorm2,
    objective = f,
    elapsed_time = Δt,
  )
  return ges
end

function my_LSR1(
  nlp::AbstractNLPModel;
  x::AbstractVector = copy(nlp.meta.x0),
  T = eltype(x),
  kwargs...,
)
  B = LSR1Operator(nlp.meta.nvar, scaling = true)::LSR1Operator{T}
  println("\n\t LSR1 trust-region method")
  return solver_TR_CG_Ab_NLP_LO_ges(nlp, B; x = x, kwargs...)
end

function my_LBFGS(
  nlp::AbstractNLPModel;
  x::AbstractVector = copy(nlp.meta.x0),
  T = eltype(x),
  kwargs...,
)
  B = LBFGSOperator(nlp.meta.nvar, scaling = true)::LBFGSOperator{T}
  println("\n\t LBFGS trust-region method")
  return solver_TR_CG_Ab_NLP_LO_ges(nlp, B; x = x, kwargs...)
end