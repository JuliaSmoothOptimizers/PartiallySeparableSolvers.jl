module Mod_ab_partitioned_data

using SparseArrays
using ReverseDiff, ForwardDiff
using PartitionedStructures, ExpressionTreeForge
using ..PartiallySeparableNLPModels.Mod_common

import NLPModels: hprod, hprod!, hess

export PartitionedData
export get_n,
  get_N,
  get_vec_elt_fun,
  get_M,
  get_vec_elt_complete_expr_tree,
  get_element_expr_tree_table,
  get_index_element_tree,
  get_vec_compiled_element_gradients
export get_x, get_v, get_s, get_pg, get_pv, get_py, get_ps, get_phv, get_pB, get_fx
export set_n!,
  set_N!,
  set_vec_elt_fun!,
  set_M!,
  set_vec_elt_complete_expr_tree!,
  set_element_expr_tree_table!,
  set_index_element_tree!,
  set_vec_compiled_element_gradients!
export set_x!,
  set_v!,
  set_s!,
  set_pg!,
  set_pv!,
  set_ps!,
  set_pg!,
  set_pv!,
  set_py!,
  set_ps!,
  set_phv!,
  set_pB!,
  set_fx!

export product_part_data_x, evaluate_obj_part_data, evaluate_grad_part_data
export product_part_data_x!,
  evaluate_obj_part_data!, evaluate_y_part_data!, evaluate_grad_part_data!
export update_nlp!

abstract type PartitionedData end

@inline get_n(part_data::PartitionedData) = part_data.n

@inline get_N(part_data::PartitionedData) = part_data.N

@inline get_vec_elt_fun(part_data::PartitionedData) = part_data.vec_elt_fun

@inline get_M(part_data::PartitionedData) = part_data.M

@inline get_vec_elt_complete_expr_tree(part_data::PartitionedData) =
  part_data.vec_elt_complete_expr_tree

@inline get_vec_elt_complete_expr_tree(part_data::PartitionedData, i::Int) =
  part_data.vec_elt_complete_expr_tree[i]

@inline get_element_expr_tree_table(part_data::PartitionedData) = part_data.element_expr_tree_table

@inline get_index_element_tree(part_data::PartitionedData) = part_data.index_element_tree

@inline get_vec_compiled_element_gradients(part_data::PartitionedData) =
  part_data.vec_compiled_element_gradients

@inline get_vec_compiled_element_gradients(part_data::PartitionedData, i::Int) =
  part_data.vec_compiled_element_gradients[i]

@inline get_x(part_data::PartitionedData) = part_data.x
@inline get_v(part_data::PartitionedData) = part_data.v
@inline get_s(part_data::PartitionedData) = part_data.s
@inline get_pg(part_data::PartitionedData) = part_data.pg
@inline get_pv(part_data::PartitionedData) = part_data.pv
@inline get_py(part_data::PartitionedData) = part_data.py
@inline get_ps(part_data::PartitionedData) = part_data.ps
@inline get_phv(part_data::PartitionedData) = part_data.phv
@inline get_pB(part_data::PartitionedData) = part_data.pB
@inline get_fx(part_data::PartitionedData) = part_data.fx

@inline set_n!(part_data::PartitionedData, n::Int) = part_data.n = n
@inline set_N!(part_data::PartitionedData, N::Int) = part_data.N = N
@inline set_vec_elt_fun!(part_data::PartitionedData, vec_elt_fun::Vector{ElementFunction}) =
  part_data.vec_elt_fun .= vec_elt_fun

@inline set_M!(part_data::PartitionedData, M::Int) = part_data.M = M

@inline set_vec_elt_complete_expr_tree!(
  part_data::PartitionedData,
  vec_elt_complete_expr_tree::Vector{G},
) where {G} = part_data.vec_elt_complete_expr_tree .= vec_elt_complete_expr_tree

@inline set_element_expr_tree_table!(
  part_data::PartitionedData,
  element_expr_tree_table::Vector{Vector{Int}},
) = part_data.element_expr_tree_table .= element_expr_tree_table

@inline set_index_element_tree!(part_data::PartitionedData, index_element_tree::Vector{Int}) =
  part_data.index_element_tree .= index_element_tree

@inline set_vec_compiled_element_gradients!(
  part_data::PartitionedData,
  vec_compiled_element_gradients::Vector{ReverseDiff.CompiledTape},
) = part_data.vec_compiled_element_gradients = vec_compiled_element_gradients

@inline set_x!(part_data::PartitionedData, x::AbstractVector{Y}) where {Y <: Number} =
  part_data.x .= x

@inline set_v!(part_data::PartitionedData, v::AbstractVector{Y}) where {Y <: Number} =
  part_data.v .= v

@inline set_s!(part_data::PartitionedData, s::AbstractVector{Y}) where {Y <: Number} =
  part_data.s .= s

@inline set_pg!(
  part_data::PartitionedData,
  pg::PartitionedStructures.Elemental_pv{Y},
) where {Y <: Number} = PartitionedStructures.epv_from_epv!(part_data.pg, pg)

@inline set_pv!(
  part_data::PartitionedData,
  pv::PartitionedStructures.Elemental_pv{Y},
) where {Y <: Number} = PartitionedStructures.epv_from_epv!(part_data.pv, pv)

@inline set_py!(
  part_data::PartitionedData,
  py::PartitionedStructures.Elemental_pv{Y},
) where {Y <: Number} = PartitionedStructures.epv_from_epv!(part_data.py, py)

@inline set_ps!(
  part_data::PartitionedData,
  ps::PartitionedStructures.Elemental_pv{Y},
) where {Y <: Number} = PartitionedStructures.epv_from_epv!(part_data.ps, ps)

@inline set_pg!(part_data::PartitionedData, x::AbstractVector{Y}) where {Y <: Number} =
  PartitionedStructures.epv_from_v!(part_data.pg, x)

@inline set_pv!(part_data::PartitionedData, v::AbstractVector{Y}) where {Y <: Number} =
  PartitionedStructures.epv_from_v!(part_data.pv, v)

@inline set_py!(part_data::PartitionedData, y::AbstractVector{Y}) where {Y <: Number} =
  PartitionedStructures.epv_from_v!(part_data.py, y)

@inline set_ps!(part_data::PartitionedData, s::AbstractVector{Y}) where {Y <: Number} =
  PartitionedStructures.epv_from_v!(part_data.ps, s)

@inline set_phv!(part_data::PartitionedData, s::AbstractVector{Y}) where {Y <: Number} =
  PartitionedStructures.epv_from_v!(part_data.phv, s)

@inline set_pB!(
  part_data::PartitionedData,
  pB::PartitionedStructures.Elemental_pm{Y},
) where {Y <: Number} = part_data.pB = pB

@inline set_fx!(part_data::PartitionedData, fx::Y) where {Y <: Number} = part_data.fx = fx

update_nlp!(part_data::PartitionedData) = @error("Should not be called")

"""
    Bx = product_part_data_x(part_data::PartitionedData, v::AbstractVector{Y}) where {Y <: Number}

Return the product between the partitioned matrix `part_data.pB` and `v`.
"""
function product_part_data_x(part_data::PartitionedData, v::AbstractVector{Y}) where {Y <: Number}
  res = similar(v)
  product_part_data_x!(res, part_data, v)
  return res
end

"""
    product_part_data_x!(res::AbstractVector{Y}, part_data::PartitionedData, v::AbstractVector{Y}) where {Y <: Number}
    product_part_data_x!(epv_res::PartitionedStructures.Elemental_pv{Y}, part_data::PartitionedData, epv::PartitionedStructures.Elemental_pv{Y}) where {Y <: Number}
    product_part_data_x!(epv_res::PartitionedStructures.Elemental_pv{Y}, pB::T, epv::PartitionedStructures.Elemental_pv{Y}) where {Y <: Number, T <: PartitionedStructures.Part_mat{Y}}

Store the product between the partitioned matrix `part_data.pB` and the vector `v` in `res`.
The computation of every element matrix-vector product require two partitioned vectors.
The first partitioned vector `epv` is created from `v` and the results are stored in the second partitioned vector `epv_res` which builds the value store by `g`.
"""
function product_part_data_x!(
  res::AbstractVector{Y},
  part_data::PartitionedData,
  v::AbstractVector{Y},
) where {Y <: Number}
  pB = get_pB(part_data)
  epv = get_ps(part_data) # a first temporary partitioned vector
  PartitionedStructures.epv_from_v!(epv, v)
  epv_res = get_py(part_data) # a second temporary partitioned vector
  product_part_data_x!(epv_res, pB, epv)
  PartitionedStructures.build_v!(epv_res)
  res .= PartitionedStructures.get_v(epv_res)
  return part_data
end

@inline product_part_data_x!(
  epv_res::PartitionedStructures.Elemental_pv{Y},
  part_data::PartitionedData,
  epv::PartitionedStructures.Elemental_pv{Y},
) where {Y <: Number} = PartitionedStructures.mul_epm_epv!(epv_res, get_pB(part_data), epv)

@inline product_part_data_x!(
  epv_res::PartitionedStructures.Elemental_pv{Y},
  pB::T,
  epv::PartitionedStructures.Elemental_pv{Y},
) where {Y <: Number, T <: PartitionedStructures.Part_mat{Y}} =
  PartitionedStructures.mul_epm_epv!(epv_res, pB, epv)

"""
    fx = evaluate_obj_part_data(part_data::PartitionedData, x::AbstractVector{Y}) where {Y <: Number}

Return the partially separable objective at `x`.
"""
function evaluate_obj_part_data(
  part_data::PartitionedData,
  x::AbstractVector{Y},
) where {Y <: Number}
  set_x!(part_data, x)
  evaluate_obj_part_data!(part_data)
  return get_fx(part_data)
end

"""
    evaluate_obj_part_data!(part_data::PartitionedData)

Compute the partially separable objective function, as a sum of element functions, at `part_data.x`.
Its objective value is stored in `part_data.fx` .
"""
function evaluate_obj_part_data!(part_data::PartitionedData)
  set_pv!(part_data, get_x(part_data))
  index_element_tree = get_index_element_tree(part_data)
  N = get_N(part_data)
  acc = 0
  for i = 1:N
    elt_expr_tree = get_vec_elt_complete_expr_tree(part_data, index_element_tree[i])
    fix = ExpressionTreeForge.evaluate_expr_tree(
      elt_expr_tree,
      PartitionedStructures.get_eev_value(get_pv(part_data), i),
    )
    acc += fix
  end
  set_fx!(part_data, acc)
  return part_data
end

"""
    evaluate_y_part_data!(part_data::PartitionedData, x::AbstractVector{Y}, s::AbstractVector{Y}) where {Y <: Number}
    evaluate_y_part_data!(part_data::PartitionedData, s::AbstractVector{Y}) where {Y <: Number}

Compute element gradients differences ∇̂fᵢ(x+s)-∇̂fᵢ(x) for each element function. 
Store the results in `part_data.py`.
When `x` is ommited, consider `part_data.x` as `x` and `part_data.pg` as the partitioned gradient at `x`.
"""
function evaluate_y_part_data!(
  part_data::PartitionedData,
  x::AbstractVector{Y},
  s::AbstractVector{Y},
) where {Y <: Number}
  set_x!(part_data, x)
  evaluate_grad_part_data!(part_data)
  evaluate_y_part_data!(part_data, s)
  return part_data
end

function evaluate_y_part_data!(part_data::PartitionedData, s::AbstractVector{Y}) where {Y <: Number}
  set_s!(part_data, s)
  set_py!(part_data, get_pg(part_data))
  PartitionedStructures.minus_epv!(get_py(part_data))
  set_x!(part_data, get_x(part_data) + s)
  evaluate_grad_part_data!(part_data)
  PartitionedStructures.add_epv!(get_pg(part_data), get_py(part_data))
  return part_data
end

"""
    gradient = evaluate_grad_part_data(part_data::PartitionedData, x::AbstractVector{Y}) where {Y <: Number}

Return the gradient vector `g` at `x`.
After the computation of the element gradients (stored in `part_data.pg`), `g` is built by accumulating their contributions.
"""
function evaluate_grad_part_data(
  part_data::PartitionedData,
  x::AbstractVector{Y},
) where {Y <: Number}
  g = similar(x)
  evaluate_grad_part_data!(g, part_data, x)
  return g
end

"""
    evaluate_grad_part_data!(g::AbstractVector{Y}, part_data::PartitionedData, x::AbstractVector{Y}) where {Y <: Number}
    evaluate_grad_part_data!(part_data::PartitionedData)

Evaluate the gradient at `x` in place.
It computes and store the element gradients in `part_data.pg` and accumulate their contributions in `g`.
When `g` and `x` are omitted, consider that `part_data.g` and `part_data.x` are `g` and `x`, respectively.
"""
function evaluate_grad_part_data!(
  g::AbstractVector{Y},
  part_data::PartitionedData,
  x::AbstractVector{Y},
) where {Y <: Number}
  x != get_x(part_data) && set_x!(part_data, x)
  evaluate_grad_part_data!(part_data)
  g .= PartitionedStructures.get_v(get_pg(part_data))
  return g
end

function evaluate_grad_part_data!(part_data::PartitionedData)
  set_pv!(part_data, get_x(part_data))
  pg = get_pg(part_data)
  index_element_tree = get_index_element_tree(part_data)
  N = get_N(part_data)
  for i = 1:N
    compiled_tape = get_vec_compiled_element_gradients(part_data, index_element_tree[i])
    Uix = PartitionedStructures.get_eev_value(get_pv(part_data), i)
    gi = PartitionedStructures.get_eev_value(get_pg(part_data), i)
    ReverseDiff.gradient!(gi, compiled_tape, Uix)
  end
  PartitionedStructures.build_v!(pg)
  return part_data
end

"""
    hv = hprod(part_data::PartitionedData, x::AbstractVector, v::AbstractVector)

Return `hv::Vector` the product between the partitioned hessian ∇²f(`x`) and the vector `v`.
"""
function hprod(part_data::PartitionedData, x::AbstractVector, v::AbstractVector)
  hv = similar(x)
  hprod!(part_data, x, v, hv)
end

"""
    hprod!(part_data::PartitionedData, x::AbstractVector, v::AbstractVector, hv::AbstractVector)

Set `hv::Vector` to the product between the partitioned hessian ∇²f(`x`) and the vector `v`.
"""
function hprod!(
  part_data::PartitionedData,
  x::AbstractVector,
  v::AbstractVector,
  hv::AbstractVector,
)
  set_ps!(part_data, x)
  set_pv!(part_data, v)

  index_element_tree = get_index_element_tree(part_data)
  N = get_N(part_data)
  ∇f(x; f) = ReverseDiff.gradient(f, x)
  ∇²fv!(x, v, hv; f) = ForwardDiff.derivative!(hv, t -> ∇f(x + t * v; f), 0)

  for i = 1:N
    complete_tree = get_vec_elt_complete_expr_tree(part_data, index_element_tree[i])
    elf_fun = ExpressionTreeForge.evaluate_expr_tree(complete_tree)

    Uix = PartitionedStructures.get_eev_value(get_ps(part_data), i)
    Uiv = PartitionedStructures.get_eev_value(get_pv(part_data), i)

    hvi = PartitionedStructures.get_eev_value(get_phv(part_data), i)
    ∇²fv!(Uix, Uiv, hvi; f = elf_fun)
  end
  PartitionedStructures.build_v!(get_phv(part_data))
  hv .= PartitionedStructures.get_v(get_phv(part_data))
end

end
