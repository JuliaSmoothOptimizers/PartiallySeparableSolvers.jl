using LinearAlgebra

# using Test, BenchmarkTools, ProfileView, InteractiveUtils


"""
    update_BFGS(Δx, y, B, B_1 )
function that builde the next approximation of the Hessian, which will be stored in B_1. The result is based from Δx, y and B; respectively the difference
between 2 points, the difference of the gradient of the associate points, and the previous approximation of the Hessian. The approximation is made according
to the SR1 update method.
"""
    function update_SR1!(Δx :: AbstractVector{Y},
                        y :: AbstractVector{Y},
                        B :: AbstractArray{Y,2}, B_1 :: AbstractArray{Y,2}) where Y <: Number
        n = length(Δx)
        ω = 1e-8
        @inbounds @fastmath v = y - B * Δx :: Vector{Y}

        @fastmath cond_left = abs( dot(Δx', v))
        @fastmath cond_right = ω * norm(Δx,2) * norm(v,2)
        @fastmath cond = (cond_left > cond_right) :: Bool

        if cond
            @fastmath num = Array{Y,2}( v * v')
            @fastmath den = dot(v', Δx) :: Y
            @fastmath num_den = num/den
            @fastmath @inbounds B_1 .= (B .+ num_den) :: Array{Y,2}
        else
            @inbounds B_1 .= B :: Array{Y,2}
        end
    end

    function update_SR1!(x :: Vector{Y}, x_1 :: Vector{Y},
                        g :: Vector{Y}, g_1 :: Vector{Y},
                        B :: Array{Y,2}, B_1 :: Array{Y,2}) where Y <: Number
        update_SR1!( x_1 - x, g_1 - g, B, B_1)
    end


"""
    update_BFGS(Δx, y, B, B_1 )
function that builde the next approximation of the Hessian, which will be stored in B_1. Rhe result is based from Δx, y and B; respectively the difference
between 2 points, the difference of the gradient of the associate points, and the previous approximation of the Hessian. The approximation is made according
to the BFGS update method.
"""
    function update_BFGS!(Δx :: AbstractVector{Y}, #difference between to points
                        y :: AbstractVector{Y}, #diffrence of the gradient between each point
                        B :: AbstractArray{Y,2}, #current approcimation of the Hessian
                        B_1 :: AbstractArray{Y,2}) where Y <: Number #Array that will store the next approximation of the Hessian

        if (dot(Δx, y) > 0 )
            @fastmath α = 1 / dot(y, Δx)
            @fastmath β = - (1 / dot(Δx, B * Δx) )
            @fastmath u = y
            @fastmath v = B * Δx
            @fastmath terme1 = (α * u * u')
            @fastmath terme2 = (β * v * v')
            @fastmath @inbounds B_1 .= (B .+ terme1 .+ terme2) :: Array{Y,2}
						return 1
        else
            @inbounds B_1 .= B :: Array{Y,2}
						return 0
        end
    end

    function update_BFGS!(x :: AbstractVector{Y}, x_1 :: AbstractVector{Y},
                        g :: AbstractVector{Y}, g_1 :: AbstractVector{Y},
                        B :: AbstractArray{Y,2}, B_1 :: AbstractArray{Y,2}) where Y <: Number
        update_BFGS!(x_1 - x, g_1 - g, B, B_1)
    end


"""
    update_SPS_SR1(sps, Bₖ, Bₖ₊₁, yₖ, sₖ)
update the Hessian approximation Bₖ using the SR1 method, according to the sps partially separable structre. To make
the update, we need the grad_vector y and the vector s. B, B_1 and y use structure linked with the partially separable structure stored in sps.
"""
    function update_SPS_SR1!(sps :: PartiallySeparableNLPModels.SPS{T},
                             B :: PartiallySeparableNLPModels.Hess_matrix{Y},
                             B_1 :: PartiallySeparableNLPModels.Hess_matrix{Y},
                             y :: PartiallySeparableNLPModels.grad_vector{Y},
                             s :: AbstractVector{Y}) where T where Y <: Number
        l_elmt_fun = length(sps.structure)
         # @Threads.threads for i in 1:l_elmt_fun
         for i in 1:l_elmt_fun
            @inbounds s_elem = Array(view(s, sps.structure[i].used_variable))
            @inbounds y_elem = y.arr[i].g_i
            @inbounds B_elem = B.arr[i].elmt_hess
            @inbounds B_elem_1 = B_1.arr[i].elmt_hess
            update_SR1!(s_elem, y_elem, B_elem, B_elem_1)
        end
    end


"""
    update_SPS_BFGS(sps, Bₖ, Bₖ₊₁, yₖ, sₖ)
update the Hessian approximation Bₖ using the SR1 method, according to the sps partially separable structre. To make
the update, we need the grad_vector y and the vector s. B, B_1 and y use structure linked with the partially separable structure stored in sps.
"""
    function update_SPS_BFGS!(sps :: PartiallySeparableNLPModels.SPS{T},
                              B :: PartiallySeparableNLPModels.Hess_matrix{Y},
                              B_1 :: PartiallySeparableNLPModels.Hess_matrix{Y},
                              y :: PartiallySeparableNLPModels.grad_vector{Y},
                              s :: AbstractVector{Y}) where T where Y <: Number
        l_elmt_fun = length(sps.structure)
         # @Threads.threads for i in 1:l_elmt_fun
				 acc = 0
         for i in 1:l_elmt_fun
            @inbounds s_elem = Array(view(s, sps.structure[i].used_variable))
            @inbounds y_elem = y.arr[i].g_i
            @inbounds B_elem = B.arr[i].elmt_hess
            @inbounds B_elem_1 = B_1.arr[i].elmt_hess
            acc += update_BFGS!(s_elem, y_elem, B_elem, B_elem_1)
        end
				println("PBFGS (old) : $(acc)/$(l_elmt_fun)")
    end



    """
        update_SPS_mix_SR1_BFGS(sps, Bₖ, Bₖ₊₁, yₖ, sₖ)
    The update is based on the grad_vector y and the vector s. B, B_1 and y use structure linked with the partially separable structure stored in sps.
    We update the Hessian approximation Bₖ according with partially separable structure. We update each Bᵢₖ,if it is convex we use
    BFGS rather than SR1.
    """
    function update_SPS_mix_SR1_BFGS!(sps :: PartiallySeparableNLPModels.SPS{T},
                                      B :: PartiallySeparableNLPModels.Hess_matrix{Y},
                                      B_1 :: PartiallySeparableNLPModels.Hess_matrix{Y},
                                      y :: PartiallySeparableNLPModels.grad_vector{Y},
                                      s :: AbstractVector{Y}) where T where Y <: Number
    l_elmt_fun = length(sps.structure)
     # @Threads.threads for i in 1:l_elmt_fun
    for i in 1:l_elmt_fun
        elmt_fun = sps.structure[i]
        status = PartiallySeparableNLPModels.get_convexity_status(elmt_fun)
        s_elem = Array(view(s, elmt_fun.used_variable))
        @inbounds y_elem = y.arr[i].g_i
        @inbounds B_elem = B.arr[i].elmt_hess
        @inbounds B_elem_1 = B_1.arr[i].elmt_hess
        if CalculusTreeTools.is_convex(status)
            update_BFGS!(s_elem, y_elem, B_elem, B_elem_1)
        else
            update_SR1!(s_elem, y_elem, B_elem, B_elem_1)
        end
    end
end
