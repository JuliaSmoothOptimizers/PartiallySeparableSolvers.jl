module Mod_TR_CG_part_data

	using CalculusTreeTools, PartitionedStructures, PartiallySeparableNLPModels
	using LinearAlgebra, LinearAlgebra.BLAS, LinearOperators, NLPModels, Krylov
	using Printf, SolverCore, SolverTools
	
	export generic_algorithm_wrapper

	function generic_algorithm_wrapper(nlp :: N, part_data :: P;
		max_eval :: Int=10000,
		max_iter::Int=10000,
		start_time::Float64=time(),
		max_time::Float64=30.0,
		ϵ::Float64= 1e-6,
		name::String="Undefined_TR_CG",
		kwargs...) where {N <: AbstractNLPModel, P <: PartiallySeparableNLPModels.PartitionedData}

		x₀ = get_x(part_data)
		n = get_n(part_data)
		∇f₀ = evaluate_grad_part_data(part_data, x₀)
		∇fNorm2 = norm(∇f₀,2)

		println("Start: "*name)
		(x,iter) = TR_CG_PD(part_data; max_eval=max_eval, max_iter=max_iter, max_time=max_time, ∇f₀=∇f₀, kwargs...)

		Δt = time() - start_time
		f = evaluate_obj_part_data(part_data, x)
		g = evaluate_grad_part_data(part_data, x)
		nrm_grad = norm(g,2)
		

		absolute(n,gₖ,ϵ) = norm(gₖ,2) < ϵ
		relative(n,gₖ,ϵ,∇fNorm2) = norm(gₖ,2) < ϵ * ∇fNorm2
		_max_iter(iter, max_iter) = iter >= max_iter
		_max_time(start_time) = (time() - start_time) >= max_time

		if absolute(n,g,ϵ) || relative(n,g,ϵ,∇fNorm2)
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
		ges = GenericExecutionStats(status, nlp, solution=x, iter=iter, dual_feas = nrm_grad, objective = f, elapsed_time = Δt)
		return ges		
	end

	function TR_CG_PD(part_data :: P;
		x::AbstractVector=copy(get_x(part_data)),
		n::Int=get_n(part_data),
		max_eval::Int=10000,
		max_iter::Int=10000,
		max_time :: Float64=30.0,
		atol::Real=√eps(eltype(x)),
		rtol::Real=√eps(eltype(x)),
		start_time::Float64=time(),
		η::Float64=1e-3,
		η₁::Float64=0.75, # > η
		Δ::Float64=1.,
		ϵ::Float64=1e-6,
		δ::Float64=1.,
		ϕ::Float64=2.,
		∇f₀::AbstractVector=PartiallySeparableNLPModels.evaluate_grad_part_data(part_data, x),
		iter_print::Int64=Int(floor(max_iter/100)),
		T=eltype(x),
		kwargs...,
		) where P <: PartiallySeparableNLPModels.PartitionedData

		iter = 0 # ≈ k
		gₖ = copy(∇f₀)
		∇fNorm2 = norm(∇f₀, 2)
		sₖ = similar(x)
				
		fₖ = evaluate_obj_part_data(part_data, x)
		
		@printf "iter temps fₖ norm(gₖ,2) Δ ρₖ\n" 

		cgtol = one(T)  # Must be ≤ 1.
		cgtol = max(rtol, min(T(0.1), 9 * cgtol / 10, sqrt(∇fNorm2)))
		
		ρₖ = -1
		B = LinearOperators.LinearOperator(T,n, n, true, true, ((res,v) -> PartiallySeparableNLPModels.product_part_data_x!(res, part_data, v)) )
		
		# stop condition
		absolute(n,gₖ,ϵ) = norm(gₖ,2) > ϵ
		relative(n,gₖ,ϵ,∇fNorm2) = norm(gₖ,2) > ϵ * ∇fNorm2
		_max_iter(iter, max_iter) = iter < max_iter
		_max_time(start_time) = (time() - start_time) < max_time
		while absolute(n,gₖ,ϵ) && relative(n,gₖ,ϵ,∇fNorm2) && _max_iter(iter, max_iter) & _max_time(start_time)# stop condition
			@printf "%3d %4g %8.1e %7.1e %7.1e" iter (time() - start_time) fₖ norm(gₖ,2) Δ
			iter += 1			
			cg_res = Krylov.cg(B, - gₖ, atol=T(atol), rtol=cgtol, radius = T(Δ), itmax=max(2*n,50))
			sₖ .= cg_res[1]  # result of the linear system solved by Krylov.cg
			
			(ρₖ, fₖ₊₁) = compute_ratio(x, fₖ, sₖ, part_data, B, gₖ) # we compute the ratio
			# step acceptance + update f,g
			if ρₖ > η
				x .= x .+ sₖ
				fₖ = fₖ₊₁
				PartiallySeparableNLPModels.update!(part_data, sₖ)
				PartiallySeparableNLPModels.evaluate_grad_part_data!(gₖ, part_data, x)
				@printf "✅\n"
			else
				fₖ = evaluate_obj_part_data(part_data, x)
				@printf "❌\n"
			end
			# trust region update
			(ρₖ >= η₁ && norm(sₖ, 2) >= 0.8*Δ) ? Δ = ϕ*Δ : Δ = Δ
			(ρₖ <= η) && (Δ = 1/ϕ*Δ)			
		end
		@printf "%3d %4g %8.1e %7.1e %7.1e %7.1e\n" iter (time() - start_time) fₖ norm(gₖ,2) Δ ρₖ

		return (x, iter)
	end

	function compute_ratio(x::AbstractVector{T}, fₖ::T, sₖ::Vector{T}, part_data::P, B::AbstractLinearOperator{T}, gₖ::AbstractVector{T}) where {T <: Number, P <: PartiallySeparableNLPModels.PartitionedData}
		mₖ₊₁ = fₖ + dot(gₖ,sₖ) + 1/2 * (dot((B*sₖ),sₖ))
		fₖ₊₁ = evaluate_obj_part_data(part_data, x+sₖ)
		ρₖ = (fₖ - fₖ₊₁)/(fₖ - mₖ₊₁)
		return (ρₖ,fₖ₊₁)
	end

end 