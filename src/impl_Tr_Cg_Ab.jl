using LinearOperators, LinearAlgebra, Krylov, NLPModels
using Printf
using SolverTools


function compute_ratio(x :: AbstractVector{Y}, f_x :: Y, s :: Vector{Y}, nlp :: AbstractNLPModel , B :: AbstractLinearOperator{Y}, g :: AbstractVector{Y})  where Y <: Number
    quad_model_s =  f_x + g' * s + 1/2 * ((B * s )' * s)  :: Y
    f_next_x = NLPModels.obj(nlp, x+s) :: Y
    num = f_x - f_next_x :: Y
    den = f_x - quad_model_s :: Y
    return (num/den, f_next_x) :: Tuple{Y,Y}
end

"""
    upgrade_TR_LO(ρₖ, xₖ , sₖ, gₖ, Bₖ, nlp, Δ)
Update the LinearOperator Bₖ using the push! if the ratio ρₖ is accepted by computing xₖ₊₁ = xₖ + sₖ and gₖ₊₁ from xₖ₊₁.
Update the radius of the Trust region Δ.
"""
function upgrade_TR_LO!( pk :: Float64, # value of the ratio
                     x_k :: AbstractVector{T}, # actual point
                     s_k :: AbstractVector{T}, # point found at the iteration k
                     g_k :: AbstractVector{T}, # array of element gradient
                     B_k :: AbstractLinearOperator{T}, # SR1 approximation
                     nlp :: AbstractNLPModel ,
                     Δ :: Float64 # radius
                     ) where T <: Number
     η = 1e-3
     η1 =  0.75
     if pk > η
         x_k .= x_k + s_k
         g_p = Vector{T}(undef,length(g_k))
         g_p .= g_k
         NLPModels.grad!(nlp, x_k, g_k)
         y_k = g_k - g_p
         push!(B_k, s_k, y_k)
     else
     end
     if pk >= η1 #now we update ∆
         if norm(s_k, 2) < 0.8 * Δ
             Δ = Δ
         else
             Δ = 2 * Δ
         end
     elseif pk <= η
         Δ = 1/2 * Δ
     end
     return Δ
end


"""
    solver_L_SR1_Ab_NLP(nlp, B, x0)
solver_L_SR1_Ab_NLP is a optimisation method using Trust region with conjugate gradient. nlp is an AbstractNLPModel, B is an
AbstractLinearOperator and x0 is an AbstractVector reprensenting the initial point of the method. The method return a tuple:
with the final point and the number of iteration performming by the method.
Even if the name suggest LSR1, B can be a BFGS Operator.
"""
function solver_TR_CG_Ab_NLP_LO(nlp :: AbstractNLPModel, B :: AbstractLinearOperator{T};
    x :: AbstractVector=copy(nlp.meta.x0),
    max_eval :: Int=10000,
    atol :: Real=√eps(eltype(x)),
    rtol :: Real=√eps(eltype(x)),
    kwargs... ) where T <: Number

    T2 = eltype(x)
    (η, Δ, ϵ, cpt) = (1e-3, 1.0, 10^-6, 0)
    n = nlp.meta.nvar

    g = Vector{T}(undef,length(x))
    ∇f₀ = Vector{T}(undef,length(x))
    NLPModels.grad!(nlp, x, ∇f₀)
    g = ∇f₀
    ∇fNorm2 = nrm2(n, ∇f₀)

    f_xk = NLPModels.obj(nlp, x)
    @printf "%3d %8.1e %7.1e %7.1e  \n" cpt f_xk norm(g,2) Δ

    cgtol = one(T2)  # Must be ≤ 1.
    cgtol = max(rtol, min(T2(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))

    while (nrm2(n,g) > ϵ) && (nrm2(n,g) > ϵ * ∇fNorm2)  && cpt < max_eval  # stop condition
        cpt = cpt + 1

        cg_res = Krylov.cg(B, - g, atol=T2(atol), rtol=cgtol, radius = Δ, itmax=max(2 * n, 50))
        sk = cg_res[1]  # result of the linear system solved by Krylov.cg

        (pk, f_temp) = compute_ratio(x, f_xk, sk, nlp, B, g) # we compute the ratio

        Δ = upgrade_TR_LO!(pk, x, sk, g, B, nlp, Δ) # we upgrade x,g,B,∆

        if  pk > η
            f_xk = f_temp
        end
        if mod(cpt,5000) == 0
            @printf "%3d %8.1e %7.1e %7.1e  \n" cpt f_xk norm(g,2) Δ
        end
    end

    if cpt < max_eval
        println("\n\n\nNous nous somme arrêté grâce à un point stationnaire !!!\n\n\n")
        println("cpt,\tf_xk,\tnorm de g,\trayon puis x en dessous ")
        @printf "%3d %8.1e %7.1e %7.1e  \n" cpt f_xk norm(g,2) Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    else
        println("\n\n\nNous nous sommes arrêté à cause du nombre d'itération max \n\n\n ")
        println("cpt,\tf_xk,\tnorm de g,\trayon puis x en dessous ")
        @printf "%3d %8.1e %7.1e %7.1e  \n" cpt f_xk norm(g,2) Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    end

    return (x, cpt)
end

#=
MISE EN FORME DES RESULTATS
=#

function solver_TR_CG_Ab_NLP_LO_res(nlp :: AbstractNLPModel, B :: AbstractLinearOperator{T}; max_eval :: Int=10000, kwargs...) where T <: Number
    Δt = @timed ((x_final,cpt) = solver_TR_CG_Ab_NLP_LO(nlp, B; max_eval=max_eval, kwargs...))
    x_init = nlp.meta.x0
    nrm_grad = norm(NLPModels.grad(nlp, x_final),2)
    nrm_grad_init = norm(NLPModels.grad(nlp, x_init),2)
    if nrm_grad < nrm_grad_init*1e-6 || nrm_grad < 1e-6
        status = :first_order
    elseif cpt >= max_eval
        status = :max_eval
    else
        status = :unknown
    end
    return GenericExecutionStats(status, nlp,
                           solution = x_final,
                           iter = cpt,  # not quite the number of iterations!
                           dual_feas = nrm_grad,
                           objective = NLPModels.obj(nlp, x_final),
                           elapsed_time = Δt[2],
                          )
end

#=
FIN DE LA PARTIE COMMMUNER
-------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------
Création des différents LinearOperators,
L-SR1 puis L-BFGS
=#

function solver_TR_CG_Ab_NLP_L_SR1(nlp :: AbstractNLPModel; x :: AbstractVector=copy(nlp.meta.x0), kwargs...)
    T = eltype(x)
    B = LSR1Operator(nlp.meta.nvar, scaling=true) :: LSR1Operator{T} #scaling=true
    println("\n\n\n\tdébut L-SR1")
    return solver_TR_CG_Ab_NLP_LO_res(nlp, B;x=x, kwargs...)
end

my_LSR1(nlp :: AbstractNLPModel;kwargs...) = solver_TR_CG_Ab_NLP_L_SR1(nlp;kwargs...)




#=
FIN DE LA PARTIE L-SR1
-------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------------------------------
Création du LinearOperators L-BFGS, et mise en forme des résultats
Une première partie sans paramètres nommés
=#



function solver_TR_CG_Ab_NLP_L_BFGS(nlp :: AbstractNLPModel; x :: AbstractVector=copy(nlp.meta.x0), kwargs...)
    T = eltype(x)
    B = LBFGSOperator(nlp.meta.nvar, scaling=true) :: LBFGSOperator{T} #scaling=true
    println("\n\n\n\tdébut L-BFGS")
    return solver_TR_CG_Ab_NLP_LO_res(nlp, B;x=x, kwargs...)
end

my_LBFGS(nlp :: AbstractNLPModel;kwargs...) = solver_TR_CG_Ab_NLP_L_BFGS(nlp;kwargs...)
