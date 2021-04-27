@enum indice::Int8 fst=1 snd=2


#Définition de la structure nécessaire pour mmon algo, c'est une deuxième version
mutable struct struct_algo{T,Y <: Number}
    #structure partiellement séparable
    sps :: PartiallySeparableNLPModel.SPS{T}
    #tuple de B
    tpl_B :: Vector{PartiallySeparableNLPModel.Hess_matrix{Y}}

    #tuple des gradients et la différence entre les gradients y
    tpl_g :: Vector{PartiallySeparableNLPModel.grad_vector{Y}}
    grad :: AbstractVector{Y}
    y :: PartiallySeparableNLPModel.grad_vector{Y}

    # tuple de x
    tpl_x :: Vector{Vector{Y}}
    #tuple des fx
    tpl_f :: Vector{Y}

    index :: indice
    #constantes pour l'algo potentiellement à supprimer plus tard
    Δ :: Float64
    η :: Float64
    η₁ :: Float64
    ϵ :: Float64

    n_eval_obj :: Int
    n_eval_grad :: Int
end


"""
    alloc_struct_algo(obj, n, type)
Alloc the structure needed for the whole Trust Region algorithm, which include the gradients vectors and Hessian approximations,
the vectors to store the points xₖ and xₖ₋₁, some constants Δ, η... and the Partially separable structure of the obj function.
"""
function alloc_struct_algo(obj :: T, n :: Int, type=Float64 :: DataType ) where T

    # détéction de la structure partiellement séparable
    sps = PartiallySeparableNLPModel.deduct_partially_separable_structure(obj,n)
    # sps = PartiallySeparableNLPModel.deduct_partially_separable_structure(obj,n) :: PartiallySeparableNLPModel.SPS{T}

    # @show sps
    # construction des structure de données nécessaire pour le gradient à l'itération k/k+1 et le différence des gradients
    construct_element_grad = (y :: PartiallySeparableNLPModel.element_function -> PartiallySeparableNLPModel.element_gradient{type}(Vector{type}(zeros(type, length(y.used_variable)) )) )
    g_k = PartiallySeparableNLPModel.grad_vector{type}( construct_element_grad.(sps.structure) )
    g_k1 = PartiallySeparableNLPModel.grad_vector{type}( construct_element_grad.(sps.structure) )
    g = Vector{PartiallySeparableNLPModel.grad_vector{type}}([g_k,g_k1])
    y = PartiallySeparableNLPModel.grad_vector{type}( construct_element_grad.(sps.structure) )
    #finally a real sized gradient
    grad = Vector{type}(undef, n)

    # constructions des structures de données nécessaires pour le Hessien ou son approximation
    construct_element_hess = ( elm_fun :: PartiallySeparableNLPModel.element_function -> PartiallySeparableNLPModel.element_hessian{type}( Array{type,2}(undef, length(elm_fun.used_variable), length(elm_fun.used_variable) )) )
    B_k = PartiallySeparableNLPModel.Hess_matrix{type}(construct_element_hess.(sps.structure))
    B_k1 = PartiallySeparableNLPModel.Hess_matrix{type}(construct_element_hess.(sps.structure))
    B = Vector{PartiallySeparableNLPModel.Hess_matrix{type}}([B_k, B_k1])

    #définition des 2 points xk et x_k1
    x_k = Vector{type}(undef, n)
    x_k1 = Vector{type}(undef, n)
    x = Vector{Vector{type}}([x_k, x_k1])

    f = Vector{type}([0,0])

    # grad_y = Vector{type}(undef, n)

    index = indice(1)

    Δ = 1.0
    η = 1e-3
    η₁ =  0.75
    ϵ = 1e-6

    n_eval_obj = 0
    n_eval_grad = 0

    # allocation de la structure de donné contenant tout ce dont nous avons besoin pour l'algorithme
    # algo_struct = struct_algo(sps, (B_k, B_k1), (g_k, g_k1), grad_k, y, grad_y, (x_k, x_k1), (type)(0), (type)(0), fst) :: struct_algo{T, type}
    # algo_struct = struct_algo(sps, B, g, grad, y, x, f, index, Δ, η, η₁, ϵ) :: struct_algo{T, type}
    # return algo_struct :: struct_algo{T, type}

    algo_struct = struct_algo(sps, B, g, grad, y, x, f, index, Δ, η, η₁, ϵ, n_eval_obj, n_eval_grad)
    return algo_struct 
end




"""
init_struct_algo(struct_algo, x )
Once the structure is allocated, we can use init_struct to init the structure struct_algo at the point x, which is the initial point
of the algorithm.
When the initialisation is done we can start the algorithme.
"""
function init_struct_algo!( s_a :: struct_algo{T,Y},
														x_k :: AbstractVector{Y}) where T where Y <: Number

    s_a.index = fst
    s_a.tpl_x[Int(s_a.index)] = x_k

    PartiallySeparableNLPModel.evaluate_SPS_gradient!(s_a.sps, s_a.tpl_x[Int(s_a.index)], s_a.tpl_g[Int(s_a.index)])
    PartiallySeparableNLPModel.build_gradient!(s_a.sps, s_a.tpl_g[Int(s_a.index)], s_a.grad)
    s_a.n_eval_grad += 1

    # PartiallySeparableNLPModel.struct_hessian!(s_a.sps, s_a.x_k, s_a.B_k)
    PartiallySeparableNLPModel.id_hessian!(s_a.sps, s_a.tpl_B[Int(s_a.index)])

    s_a.tpl_f[Int(s_a.index)] = PartiallySeparableNLPModel.evaluate_SPS(s_a.sps, s_a.tpl_x[Int(s_a.index)])
    s_a.n_eval_obj += 1

end

"""
    other_index(i)
 i ∈ {1,2}. We use a trick to avoid useless copy of gradients and Hessians in the structure of the algorithm. So we need and index
 in the structure of the algorithm. The function other_index return the other index.
 The index is an enumerate type.
 If the index is currently to 1 (fst) in the structure s_a then other_index(s_a) return 2 (snd), and in the other case return 1  (fst)
"""
function other_index( s_a :: struct_algo{T,Y}) where T where Y <: Number
    if s_a.index == fst
        return snd :: indice
    else
        return fst :: indice
    end
end


"""
    approx_quad(struct_algo, pₖ   )
return the quadratic approximation m(pₖ) = fₖ + gₖᵀpₖ + 1/2.pₖᵀBpₖ. The values of  fₖ, gₖ and Bₖ are stored inside struct_algo.
"""
function approx_quad(s_a :: struct_algo{T,Y}, x :: AbstractVector{Y}) where T where Y <: Number
    s_a.tpl_f[Int(s_a.index)] + s_a.grad' * x  +  1/2 *  PartiallySeparableNLPModel.product_matrix_sps(s_a.sps, s_a.tpl_B[Int(s_a.index)], x)' * x
end

"""
    compute_ratio(struct_algo, sₖ )
Compute the ratio :   (fₖ - fₖ₊₁)/(mₖ(0)-mₖ(sₖ)) , all the data about fₖ, mₖ is stored in struct_algo
"""
function compute_ratio(s_a :: struct_algo{T,Y}, s_k :: AbstractVector{Y}) where T where Y <: Number
    fxₖ = s_a.tpl_f[Int(s_a.index)] :: Y
    fxₖ₊₁ = PartiallySeparableNLPModel.evaluate_SPS(s_a.sps, s_a.tpl_x[Int(s_a.index)] + s_k) :: Y
    s_a.n_eval_obj += 1
    quadratic_approximation = approx_quad(s_a, s_k) :: Y
    num = fxₖ - fxₖ₊₁ :: Y
    den = fxₖ - quadratic_approximation :: Y
    ρₖ = num/den :: Y
    return (ρₖ, fxₖ₊₁) :: Tuple{Y,Y}
end


#=
FIN DE LA PARTIE COMMUNE
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
P-SR1
=#

"""
    update_PSR1!(struct_algo, B)
This function perform a step of a Trust-Region method using a conjuguate-gragient method to solve the sub-problem of the Trust-Region.
B is the LinearOperator needed by the cg (conjuguate-gragient method). struct_algo stored all the data relative to the problem and is modified if step is taken .
"""
function update_PSR1!(s_a :: struct_algo{T,Y}, B :: LinearOperator{Y};
    atol :: Real=√eps(s_a.tpl_x[Int(s_a.index)]),
    rtol :: Real=√eps(s_a.tpl_x[Int(s_a.index)]),
    kwawrgs...
    ) where T where Y <: Number

    n = s_a.sps.n_var
    (s_k, info) = Krylov.cg(B, - s_a.grad, atol=atol, rtol=rtol, radius = s_a.Δ, itmax=max(2 * n, 50)) :: Tuple{Array{Y,1},Krylov.SimpleStats{Y}}

    (ρₖ, fxₖ₊₁) = compute_ratio(s_a, s_k :: Vector{Y})
    if ρₖ > s_a.η #= on accepte le nouveau point =#
        s_a.index = other_index(s_a)
        s_a.tpl_f[Int(s_a.index)] = fxₖ₊₁
        s_a.tpl_x[Int(s_a.index)] = s_a.tpl_x[Int(other_index(s_a))] + s_k

        PartiallySeparableNLPModel.evaluate_SPS_gradient!(s_a.sps, s_a.tpl_x[Int(s_a.index)], s_a.tpl_g[Int(s_a.index)])
        PartiallySeparableNLPModel.build_gradient!(s_a.sps, s_a.tpl_g[Int(s_a.index)], s_a.grad)
        PartiallySeparableNLPModel.minus_grad_vec!(s_a.tpl_g[Int(s_a.index)], s_a.tpl_g[Int(other_index(s_a))], s_a.y)
        s_a.n_eval_grad += 1

        update_SPS_SR1!(s_a.sps, s_a.tpl_B[Int(other_index(s_a))], s_a.tpl_B[Int(s_a.index)], s_a.y, s_k) #on obtient notre nouveau B_k
    else #= println("changement par référence, la structure ne bouge donc pas") =#
    end

    if ρₖ >= s_a.η₁
        if norm(s_k) < 0.8 * s_a.Δ #= le rayon ne bouge pas les conditions sur la norme ne sont pas satisfaites =#
            s_a.Δ = s_a.Δ
        else   #= le rayon augmenter=#
            s_a.Δ = s_a.Δ * 2
        end
    elseif ρₖ <= s_a.η  #=le rayon diminue=#
        s_a.Δ = 1/2 * s_a.Δ
    else  #= cas ou nous faisons ni une bonne ni une mauvais approximation=#
        s_a.Δ = s_a.Δ
    end
end


#fonction traitant le coeur de l'algorithme, réalise principalement la boucle qui incrémente un compteur et met à jour la structure d'algo par effet de bord
# De plus on effectue tous les affichage par itération dans cette fonction raison des printf
iterations_TR_PSR1!(s_a :: struct_algo{T,Y}, vec_cpt :: Vector{Int}, start_time :: Float64; kwargs...) where T where Y <: Number = begin  (cpt, elapsed_time) = iterations_TR_PSR1!(s_a, start_time; kwargs...); vec_cpt[1] = cpt ; return (vec_cpt[1], elapsed_time) end
function iterations_TR_PSR1!(s_a :: struct_algo{T,Y},
    start_time :: Float64;
    max_eval :: Int=10000,
    atol :: Real=√eps(eltype(s_a.tpl_x[Int(s_a.index)])),
    rtol :: Real=√eps(eltype(s_a.tpl_x[Int(s_a.index)])),
    max_time :: Float64=30.0,
    kwargs... )  where T where Y <: Number

    cpt = 1 :: Int64
    n = s_a.sps.n_var
    elapsed_time = 0.0

    T2 = eltype(s_a.tpl_x[Int(s_a.index)])
    ∇f₀ = s_a.grad
    ∇fNorm2 = nrm2(n, ∇f₀)
    cgtol = one(T2)  # Must be ≤ 1.
    cgtol = max(rtol, min(T2(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))

    opB(s :: struct_algo{T,Y}) = LinearOperators.LinearOperator(n, n, true, true, x -> PartiallySeparableNLPModel.product_matrix_sps(s.sps, s.tpl_B[Int(s.index)], x) ) :: LinearOperator{Y}
    @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n" cpt s_a.tpl_f[Int(s_a.index)] norm(s_a.grad,2) s_a.Δ

    while ( (norm(s_a.grad,2) > s_a.ϵ ) && (norm(s_a.grad,2) > s_a.ϵ * ∇fNorm2)  &&  s_a.n_eval_obj < max_eval ) && elapsed_time < max_time
        update_PSR1!(s_a, opB(s_a); atol=atol, rtol=cgtol)
        cpt = cpt + 1
        if mod(cpt,500) == 0
            @printf "\n%3d \t%8.1e \t%7.1e \t%7.1e \n" cpt s_a.tpl_f[Int(s_a.index)] norm(s_a.grad,2) s_a.Δ
        end
        elapsed_time = time() - start_time
    end

    if cpt < max_eval
        println("\n\n\nNous nous somme arrêté grâce à un point stationnaire PSR1 !!!")
        println("cpt,\tf_xk,\tnorm de g,\trayon puis x en dessous ")
        @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n\n\n" cpt s_a.tpl_f[Int(s_a.index)]  norm(s_a.grad,2) s_a.Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    else
        println("\n\n\nNous nous sommes arrêté à cause du nombre d'itération max PSR1")
        println("cpt,\tf_xk,\tnorm de g,\trayon ")
        @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n\n\n" cpt s_a.tpl_f[Int(s_a.index)]  norm(s_a.grad,2) s_a.Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    end

    return (cpt, elapsed_time)
end



"""
    solver_TR_PSR1!(p)
Trust region method using the gradient conjugate method and the Partially Separable Structure of the problem of the parameters. This method use
SR1 approximation.
"""

solver_TR_PSR1!(m :: T; kwargs...) where T <: AbstractNLPModel = _solver_TR_PSR1!(m; kwargs... )
solver_TR_PSR1!(obj_Expr :: T, n :: Int, x_init :: AbstractVector{Y}, type=Float64 :: DataType; kwargs... ) where T where Y <: Number = _solver_TR_PSR1!(obj_Expr, n, x_init, type; kwargs...)
function _solver_TR_PSR1!(obj_Expr :: T, n :: Int, x_init :: AbstractVector{Y}, type=Float64 :: DataType; kwargs... ) where T where Y <: Number

    s_a = alloc_struct_algo(obj_Expr, n :: Int, type :: DataType )
    init_struct_algo!(s_a, x_init)
    pointer_cpt = Array{Int}([0])
    start_time = time()


    # try
        (cpt, elapsed_time) = iterations_TR_PSR1!(s_a, pointer_cpt, start_time; kwargs...)
    # catch e
    #     @show e
    #     return (zeros(n),s_a,0)
    # end

    cpt = pointer_cpt[1]
    x_final = s_a.tpl_x[Int(s_a.index)]
    return (x_final, s_a, cpt, elapsed_time) :: Tuple{Vector{Y}, struct_algo, Int, Float64}
end

function _solver_TR_PSR1_2!(m :: Z, obj_Expr :: T, n :: Int, type:: DataType, x_init :: AbstractVector{Y};
    max_eval :: Int=10000,
    max_time :: Float64=30.0,
    kwargs...) where T where Y where Z <: AbstractNLPModel
    (x_final, s_a, cpt, elapsed_time) = solver_TR_PSR1!(obj_Expr, n, x_init, type; max_eval=max_eval, max_time=max_time, kwargs...)
    nrm_grad = norm(NLPModels.grad(m, x_final),2)
    nrm_grad_init = norm(NLPModels.grad(m, x_init),2)
    if nrm_grad < nrm_grad_init*1e-6 || nrm_grad < 1e-6
        status = :first_order
    elseif cpt >= max_eval
        status = :max_eval
    elseif elapsed_time > max_time
        status = :max_time
    else
        status = :unknown
    end
    m.counters.neval_obj = s_a.n_eval_obj
    m.counters.neval_grad = s_a.n_eval_grad
		# attention renvoie de s_a ajouter pour pouvoir former les matrices.
    return s_a, GenericExecutionStats(status, m,
                           solution = x_final,
                           iter = cpt,  # not quite the number of iterations!
                           dual_feas = nrm_grad,
                           objective = NLPModels.obj(m, x_final),
                           elapsed_time = elapsed_time,
                          )
end


function _solver_TR_PSR1!( model_JUMP :: T; x :: AbstractVector=copy(model_JUMP.meta.x0), kwargs...) where T <: AbstractNLPModel 
    model = model_JUMP.eval.m
    evaluator = JuMP.NLPEvaluator(model)
    MathOptInterface.initialize(evaluator, [:ExprGraph])
    obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
    n = model.moi_backend.model_cache.model.num_variables_created
    T2 = eltype(x)
    _solver_TR_PSR1_2!(model_JUMP, obj_Expr, n, T2, x; kwargs...)
end


function _solver_TR_PSR1!( adnlp :: RADNLPModel; x0 :: AbstractVector=copy(adnlp.meta.x0), kwargs...)
	n = length(x0)
	ModelingToolkit.@variables x[1:n]
	fun = adnlp.f(x)
	ex = CalculusTreeTools.transform_to_expr_tree(fun)
	T2 = eltype(x0)
	_solver_TR_PSR1_2!(adnlp, ex, n, T2, x0; kwargs...)
end



#=
FIN DE P-SR1
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
P-BFGS
=#




"""
    update_PBGS!(struct_algo, B)
This function perform a step of a Trust-Region method using a conjuguate-gragient method to solve the sub-problem of the Trust-Region.
B is the LinearOperator needed by the cg (conjuguate-gragient method). struct_algo stored all the data relative to the problem and is modified if step is taken .
"""
function update_PBGS!(s_a :: struct_algo{T,Y}, B :: LinearOperator{Y};
    atol :: Real=√eps(s_a.tpl_x[Int(s_a.index)]),
    rtol :: Real=√eps(s_a.tpl_x[Int(s_a.index)]),
    kwawrgs...
    ) where T where Y <: Number

    n = s_a.sps.n_var
    (s_k, info) = Krylov.cg(B, - s_a.grad,atol=atol, rtol=rtol, radius = s_a.Δ, itmax=max(2 * n, 50)) :: Tuple{Array{Y,1},Krylov.SimpleStats{Y}}
    (ρₖ, fxₖ₊₁) = compute_ratio(s_a, s_k :: Vector{Y})
    if ρₖ > s_a.η #= on accepte le nouveau point =#
        s_a.index = other_index(s_a)
        s_a.tpl_f[Int(s_a.index)] = fxₖ₊₁
        s_a.tpl_x[Int(s_a.index)] = s_a.tpl_x[Int(other_index(s_a))] + s_k

        PartiallySeparableNLPModel.evaluate_SPS_gradient!(s_a.sps, s_a.tpl_x[Int(s_a.index)], s_a.tpl_g[Int(s_a.index)])
        PartiallySeparableNLPModel.build_gradient!(s_a.sps, s_a.tpl_g[Int(s_a.index)], s_a.grad)
        PartiallySeparableNLPModel.minus_grad_vec!(s_a.tpl_g[Int(s_a.index)], s_a.tpl_g[Int(other_index(s_a))], s_a.y)
        s_a.n_eval_grad += 1

        update_SPS_BFGS!(s_a.sps, s_a.tpl_B[Int(other_index(s_a))], s_a.tpl_B[Int(s_a.index)], s_a.y, s_k) #on obtient notre nouveau B_k
    else #= println("changement par référence, la structure ne bouge donc pas") =#
    end

    if ρₖ >= s_a.η₁
        if norm(s_k) < 0.8 * s_a.Δ #= le rayon ne bouge pas les conditions sur la norme ne sont pas satisfaites =#
            s_a.Δ = s_a.Δ
        else   #= le rayon augmenter=#
            s_a.Δ = s_a.Δ * 2
        end
    elseif ρₖ <= s_a.η  #=le rayon diminue=#
        s_a.Δ = 1/2 * s_a.Δ
    else  #= cas ou nous faisons ni une bonne ni une mauvais approximation=#
        s_a.Δ = s_a.Δ
    end
end


#fonction traitant le coeur de l'algorithme, réalise principalement la boucle qui incrémente un compteur et met à jour la structure d'algo par effet de bord
# De plus on effectue tous les affichage par itération dans cette fonction raison des printf
iterations_TR_PBGFS!(s_a :: struct_algo{T,Y}, vec_cpt :: Vector{Int}, start_time :: Float64; kwargs...) where T where Y <: Number = begin  (cpt, elapsed_time) = iterations_TR_PBGFS!(s_a, start_time; kwargs...); vec_cpt[1] = cpt ; return (vec_cpt[1], elapsed_time) end
function iterations_TR_PBGFS!(s_a :: struct_algo{T,Y},
        start_time :: Float64;
        max_eval :: Int=10000,
        max_time :: Float64=30.0,
        atol :: Real=√eps(eltype(s_a.tpl_x[Int(s_a.index)])),
        rtol :: Real=√eps(eltype(s_a.tpl_x[Int(s_a.index)])),
        kwargs... )  where T where Y <: Number


    cpt = 1 :: Int64
    n = s_a.sps.n_var
    elasped_time = 0.0

    T2 = eltype(s_a.tpl_x[Int(s_a.index)])
    ∇f₀ = s_a.grad
    ∇fNorm2 = nrm2(n, ∇f₀)
    cgtol = one(T2)  # Must be ≤ 1.
    cgtol = max(rtol, min(T2(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))

    opB(s :: struct_algo{T,Y}) = LinearOperators.LinearOperator(n, n, true, true, x -> PartiallySeparableNLPModel.product_matrix_sps(s.sps, s.tpl_B[Int(s.index)], x) ) :: LinearOperator{Y}
    @printf "\n%3d \t%8.1e \t%7.1e \t%7.1e \n" cpt s_a.tpl_f[Int(s_a.index)] ∇fNorm2 s_a.Δ

    while ( (norm(s_a.grad,2) > s_a.ϵ ) && (norm(s_a.grad,2) > s_a.ϵ * ∇fNorm2)  &&  s_a.n_eval_obj < max_eval ) && elasped_time < max_time
        update_PBGS!(s_a, opB(s_a); atol=atol, rtol=cgtol)
        cpt = cpt + 1
        if mod(cpt,500) == 0
            @printf "\n%3d \t%8.1e \t%7.1e \t%7.1e \n" cpt s_a.tpl_f[Int(s_a.index)] norm(s_a.grad,2) s_a.Δ
        end
        elasped_time = time() - start_time
    end

    if cpt < max_eval
        println("\n\n\nNous nous somme arrêté grâce à un point stationnaire PBFGS !!!")
        println("cpt,\tf_xk,\tnorm de g,\trayon puis x en dessous ")
        @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n\n\n" cpt s_a.tpl_f[Int(s_a.index)]  norm(s_a.grad,2) s_a.Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    else
        println("\n\n\nNous nous sommes arrêté à cause du nombre d'itération max PBFGS")
        println("cpt,\tf_xk,\tnorm de g,\trayon ")
        @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n\n\n" cpt s_a.tpl_f[Int(s_a.index)]  norm(s_a.grad,2) s_a.Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    end

    return (cpt, elasped_time)
end


"""
    solver_TR_PSR1!(model)
Trust region method using the gradient conjugate method and the Partially Separable Structure of the model from the parameters. This method use
BFGS approximation.
"""
solver_TR_PBFGS!(m :: T;  kwargs... ) where T <: AbstractNLPModel = _solver_TR_PBFGS!(m; kwargs... )
solver_TR_PBFGS!(obj_Expr :: T, n :: Int, x_init :: AbstractVector{Y}, type=Float64 :: DataType ; kwargs...) where T where Y <: Number = _solver_TR_PBFGS!(obj_Expr, n, x_init, type; kwargs...)
function _solver_TR_PBFGS!(obj_Expr :: T, n :: Int, x_init :: AbstractVector{Y}, type=Float64 :: DataType; kwargs... ) where T where Y <: Number
    s_a = alloc_struct_algo(obj_Expr, n :: Int, type :: DataType )
		@show typeof(s_a)
    init_struct_algo!(s_a, x_init)
    pointer_cpt = Array{Int}([0])
    start_time = time()
    # try
        (cpt, elapsed_time) = iterations_TR_PBGFS!(s_a, pointer_cpt, start_time; kwargs...)
    # catch e
    #     @show e
    #     return (zeros(n),s_a,-1)
    # end
    cpt = pointer_cpt[1]
    x_final = s_a.tpl_x[Int(s_a.index)]
    return (x_final, s_a, cpt, elapsed_time) :: Tuple{Vector{Y}, struct_algo, Int, Float64}
end


function _solver_TR_PBFGS_2!(m :: Z, obj_Expr :: T, n :: Int, type:: DataType, x_init :: AbstractVector{Y};
        max_eval :: Int=10000,
        max_time :: Float64=30.0,
        kwargs...) where T where Y where Z <: AbstractNLPModel

    (x_final, s_a, cpt, elapsed_time) = solver_TR_PBFGS!(obj_Expr, n, x_init, type; max_eval=max_eval, max_time = max_time, kwargs...)
    nrm_grad = norm(NLPModels.grad(m, x_final),2)
    nrm_grad_init = norm(NLPModels.grad(m, x_init),2)
    if nrm_grad < nrm_grad_init*1e-6 || nrm_grad < 1e-6
        status = :first_order
    elseif cpt >= max_eval
        status = :max_eval
    elseif elapsed_time > max_time
        status = :max_time
    else
        status = :unknown
    end
    m.counters.neval_obj = s_a.n_eval_obj
    m.counters.neval_grad = s_a.n_eval_grad
    # return GenericExecutionStats(status, m,
    return s_a, GenericExecutionStats(status, m,
                           solution = x_final,
                           iter = cpt,  # not quite the number of iterations!
                           dual_feas = nrm_grad,
                           objective = NLPModels.obj(m, x_final),
                           elapsed_time = elapsed_time,
                          )
end


function _solver_TR_PBFGS!( model_JUMP :: T; x :: AbstractVector=copy(model_JUMP.meta.x0), kwargs...) where T <: AbstractNLPModel where Y <: Number
    model = model_JUMP.eval.m
    evaluator = JuMP.NLPEvaluator(model)
    MathOptInterface.initialize(evaluator, [:ExprGraph])
    obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
    n = model.moi_backend.model_cache.model.num_variables_created
    T2 = eltype(x)
    _solver_TR_PBFGS_2!(model_JUMP, obj_Expr, n, T2, x; kwargs...)
end


function _solver_TR_PBFGS!( adnlp :: RADNLPModel; x0 :: AbstractVector=copy(adnlp.meta.x0), kwargs...) 
	n = length(x0)
	ModelingToolkit.@variables x[1:n]
	fun = adnlp.f(x)
	ex = CalculusTreeTools.transform_to_expr_tree(fun)
	CalculusTreeTools.print_tree(ex)
	T2 = eltype(x0)
	_solver_TR_PBFGS_2!(adnlp, ex, n, T2, x0; kwargs...)
end



#=
FIN DE P-BFGS
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------
P-BS
=#


"""
    update_PBS!(struct_algo, B)
This function perform a step of a Trust-Region method using a conjuguate-gragient method to solve the sub-problem of the Trust-Region.
B is the LinearOperator needed by the cg (conjuguate-gragient method). struct_algo stored all the data relative to the problem and is modified if step is taken .
"""
function update_PBS!(s_a :: struct_algo{T,Y}, B :: LinearOperator{Y};
    atol :: Real=√eps(s_a.tpl_x[Int(s_a.index)]),
    rtol :: Real=√eps(s_a.tpl_x[Int(s_a.index)]),
    kwawrgs...
    ) where T where Y <: Number

    n = s_a.sps.n_var
    (s_k, info) = Krylov.cg(B, - s_a.grad,atol=atol, rtol=rtol, radius = s_a.Δ, itmax=max(2 * n, 50)) :: Tuple{Array{Y,1},Krylov.SimpleStats{Y}}
    (ρₖ, fxₖ₊₁) = compute_ratio(s_a, s_k :: Vector{Y})
    if ρₖ > s_a.η #= on accepte le nouveau point =#
        s_a.index = other_index(s_a)
        s_a.tpl_f[Int(s_a.index)] = fxₖ₊₁
        s_a.tpl_x[Int(s_a.index)] = s_a.tpl_x[Int(other_index(s_a))] + s_k

        PartiallySeparableNLPModel.evaluate_SPS_gradient!(s_a.sps, s_a.tpl_x[Int(s_a.index)], s_a.tpl_g[Int(s_a.index)])
        PartiallySeparableNLPModel.build_gradient!(s_a.sps, s_a.tpl_g[Int(s_a.index)], s_a.grad)
        PartiallySeparableNLPModel.minus_grad_vec!(s_a.tpl_g[Int(s_a.index)], s_a.tpl_g[Int(other_index(s_a))], s_a.y)
        s_a.n_eval_grad += 1

        update_SPS_mix_SR1_BFGS!(s_a.sps, s_a.tpl_B[Int(other_index(s_a))], s_a.tpl_B[Int(s_a.index)], s_a.y, s_k) #on obtient notre nouveau B_k
    else #= println("changement par référence, la structure ne bouge donc pas") =#
    end

    if ρₖ >= s_a.η₁
        if norm(s_k) < 0.8 * s_a.Δ #= le rayon ne bouge pas les conditions sur la norme ne sont pas satisfaites =#
            s_a.Δ = s_a.Δ
        else   #= le rayon augmenter=#
            s_a.Δ = s_a.Δ * 2
        end
    elseif ρₖ <= s_a.η  #=le rayon diminue=#
        s_a.Δ = 1/2 * s_a.Δ
    else  #= cas ou nous faisons ni une bonne ni une mauvais approximation=#
        s_a.Δ = s_a.Δ
    end
end


#fonction traitant le coeur de l'algorithme, réalise principalement la boucle qui incrémente un compteur et met à jour la structure d'algo par effet de bord
# De plus on effectue tous les affichage par itération dans cette fonction raison des printf
iterations_TR_PBS!(s_a :: struct_algo{T,Y}, vec_cpt :: Vector{Int}, start_time :: Float64; kwargs...) where T where Y <: Number = begin  (cpt, elapsed_time) = iterations_TR_PBS!(s_a, start_time; kwargs...); vec_cpt[1] = cpt ; return (vec_cpt[1], elapsed_time) end
function iterations_TR_PBS!(s_a :: struct_algo{T,Y},
        start_time :: Float64;
        max_eval :: Int=10000,
        max_time :: Float64=30.0,
        atol :: Real=√eps(eltype(s_a.tpl_x[Int(s_a.index)])),
        rtol :: Real=√eps(eltype(s_a.tpl_x[Int(s_a.index)])),
        kwargs... )  where T where Y <: Number

    elapsed_time = 0.0
    cpt = 1 :: Int64
    n = s_a.sps.n_var

    T2 = eltype(s_a.tpl_x[Int(s_a.index)])
    ∇f₀ = s_a.grad
    ∇fNorm2 = nrm2(n, ∇f₀)
    cgtol = one(T2)  # Must be ≤ 1.
    cgtol = max(rtol, min(T2(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))

    opB(s :: struct_algo{T,Y}) = LinearOperators.LinearOperator(n, n, true, true, x -> PartiallySeparableNLPModel.product_matrix_sps(s.sps, s.tpl_B[Int(s.index)], x) ) :: LinearOperator{Y}
    @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n" cpt s_a.tpl_f[Int(s_a.index)] ∇fNorm2 s_a.Δ

    while ( (norm(s_a.grad,2) > s_a.ϵ ) && (norm(s_a.grad,2) > s_a.ϵ * ∇fNorm2)  &&  s_a.n_eval_obj < max_eval ) && elapsed_time < max_time
        update_PBS!(s_a, opB(s_a); atol=atol, rtol=cgtol)
        cpt = cpt + 1
        if mod(cpt,500) == 0
            @printf "\n%3d \t%8.1e \t%7.1e \t%7.1e \n" cpt s_a.tpl_f[Int(s_a.index)] norm(s_a.grad,2) s_a.Δ
        end
        elapsed_time = time() - start_time
    end

    if cpt < max_eval
        println("\n\n\nNous nous somme arrêté grâce à un point stationnaire PBS !!!")
        println("cpt,\tf_xk,\tnorm de g,\trayon puis x en dessous ")
        @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n\n\n" cpt s_a.tpl_f[Int(s_a.index)]  norm(s_a.grad,2) s_a.Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    else
        println("\n\n\nNous nous sommes arrêté à cause du nombre d'itération max PBS")
        println("cpt,\tf_xk,\tnorm de g,\trayon ")
        @printf "%3d \t%8.1e \t%7.1e \t%7.1e \n\n\n" cpt s_a.tpl_f[Int(s_a.index)]  norm(s_a.grad,2) s_a.Δ
        println("---------------------------------------------------------------------------------------------------------\n\n\n")
    end

    return (cpt, elapsed_time)
end


solver_TR_PBS!(m :: T;  kwargs... ) where T <: AbstractNLPModel = _solver_TR_PBS!(m; kwargs... )
solver_TR_PBS!(obj_Expr :: T, n :: Int, x_init :: AbstractVector{Y}, type=Float64 :: DataType ; kwargs...) where T where Y <: Number = _solver_TR_PBS!(obj_Expr, n, x_init, type; kwargs...)
function _solver_TR_PBS!(obj_Expr :: T, n :: Int, x_init :: AbstractVector{Y}, type=Float64 :: DataType; kwargs... ) where T where Y <: Number
    s_a = alloc_struct_algo(obj_Expr, n :: Int, type :: DataType )
    init_struct_algo!(s_a, x_init)
    pointer_cpt = Array{Int}([0])

    start_time = time()
    # try
        (cpt, elapsed_time) = iterations_TR_PBS!(s_a, pointer_cpt, start_time; kwargs...)
    # catch e
    #     @show e
    #     return (zeros(n),s_a,-1)
    # end
    cpt = pointer_cpt[1]
    x_final = s_a.tpl_x[Int(s_a.index)]
    return (x_final, s_a, cpt, elapsed_time) :: Tuple{Vector{Y}, struct_algo, Int, Float64}
end


function _solver_TR_PBS_2!(m :: Z, obj_Expr :: T, n :: Int, type:: DataType, x_init :: AbstractVector{Y};
        max_eval :: Int=10000,
        max_time :: Float64=30.0,
        kwargs...) where T where Y where Z <: AbstractNLPModel

    (x_final, s_a, cpt, elapsed_time) = solver_TR_PBS!(obj_Expr, n, x_init, type; max_eval=max_eval, max_time=max_time, kwargs...)
    nrm_grad = norm(NLPModels.grad(m, x_final),2)
    nrm_grad_init = norm(NLPModels.grad(m, x_init),2)
    if nrm_grad < nrm_grad_init*1e-6 || nrm_grad < 1e-6
        status = :first_order
    elseif cpt >= max_eval
        status = :max_eval
    elseif elapsed_time > max_time
        status = :max_time
    else
        status = :unknown
    end
    m.counters.neval_obj = s_a.n_eval_obj
    m.counters.neval_grad = s_a.n_eval_grad
    return s_a, GenericExecutionStats(status, m,
                           solution = x_final,
                           iter = cpt,  # not quite the number of iterations!
                           dual_feas = nrm_grad,
                           objective = NLPModels.obj(m, x_final),
                           elapsed_time = elapsed_time,
                          )
end


function _solver_TR_PBS!( model_JUMP :: T; x :: AbstractVector=copy(model_JUMP.meta.x0), kwargs...) where T <: AbstractNLPModel where Y <: Number
    model = model_JUMP.eval.m
    evaluator = JuMP.NLPEvaluator(model)
    MathOptInterface.initialize(evaluator, [:ExprGraph])
    obj_Expr = MathOptInterface.objective_expr(evaluator) :: Expr
    n = model.moi_backend.model_cache.model.num_variables_created
    T2 = eltype(x)
    _solver_TR_PBS_2!(model_JUMP, obj_Expr, n, T2, x; kwargs...)
end


function _solver_TR_PBS!( adnlp :: RADNLPModel; x0 :: AbstractVector=copy(adnlp.meta.x0), kwargs...)
	n = length(x0)
	ModelingToolkit.@variables x[1:n]
	fun = adnlp.f(x)
	ex = CalculusTreeTools.transform_to_expr_tree(fun)
	T2 = eltype(x0)
	_solver_TR_PBS_2!(adnlp, ex, n, T2, x0; kwargs...)
end
