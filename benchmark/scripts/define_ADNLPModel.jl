using ADNLPModels

start_ones(n::Int) = ones(n)

function arwhead(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("arwhead: number of variables must be ≥ 2")
  n = max(2, n)

  return sum((x[i]^2 + x[n]^2)^2 - 4 * x[i] + 3 for i = 1:(n - 1))
end
start_arwhead(n::Int) = ones(n)
arwhead_ADNLPModel(n::Int = 100) =
  ADNLPModel(arwhead, start_arwhead(n), name = "arwhead " * string(n) * " variables")

function bdqrtic(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 5 && @warn("bdqrtic: number of variables must be ≥ 5")
  n = max(5, n)

  return sum(
    (3 - 4 * x[i])^2 + (x[i]^2 + 2 * x[i + 1]^2 + 3 * x[i + 2]^2 + 4 * x[i + 3]^2 + 5 * x[n]^2)^2
    for i = 1:(n - 4)
  )
end
start_bdqrtic(n::Int) = ones(n)
bdqrtic_ADNLPModel(n::Int = 100) where {Y <: Number} =
  ADNLPModel(bdqrtic, start_bdqrtic(n), name = "bdqrtic " * string(n) * " variables")

function brybnd(x::AbstractVector{Y}; ml::Int = 5, mu::Int = 1) where {Y <: Number}
  n = length(x)

  sum(
    (
      x[i] * (2 + 5 * x[i]^2) + 1 -
      sum(x[j] * (1 + x[j]) for j = max(1, i - ml):min(n, i + mu) if j != i)
    )^2 for i = 1:n
  )
end
start_brybnd(n::Int) = (x -> -1 * x).(ones(n))
brybnd_ADNLPModel(n::Int = 100) =
  ADNLPModel(brybnd, start_bdqrtic(n), name = "brybnd " * string(n) * " variables")

function chainwoo(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  (n % 4 == 0) || @warn("chainwoo: number of variables adjusted to be a multiple of 4")
  n = 4 * max(1, div(n, 4))

  1.0 + sum(
    100 * (x[2 * i] - x[2 * i - 1]^2)^2 +
    (1 - x[2 * i - 1])^2 +
    90 * (x[2 * i + 2] - x[2 * i + 1]^2)^2 +
    (1 - x[2 * i + 1])^2 +
    10 * (x[2 * i] + x[2 * i + 2] - 2)^2 +
    0.1 * (x[2 * i] - x[2 * i + 2])^2 for i = 1:(div(n, 2) - 1)
  )
end
function start_chainwoo(n::Int)
  x0 = (x -> -2 * x).(ones(n))
  x0[1] = -3
  x0[2] = -1
  x0[3] = -3
  x0[4] = -1
  return x0
end
chainwoo_ADNLPModel(n::Int = 100) =
  ADNLPModel(chainwoo, start_chainwoo(n), name = "chainwoo " * string(n) * " variables")

function cosine(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("cosine: number of variables must be ≥ 2")
  n = max(2, n)
  sum(cos(x[i]^2 - 0.5 * x[i + 1]) for i = 1:(n - 1))
end
start_cosine(n::Int) = start_ones(n)
cosine_ADNLPModel(n::Int = 100) =
  ADNLPModel(cosine, start_cosine(n), name = "cosine " * string(n) * " variables")

function cragglvy(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("cragglvy: number of variables must be ≥ 2")
  n = max(2, n)

  sum(
    (exp(x[2 * i - 1]) - x[2 * i])^4 +
    100 * (x[2 * i] - x[2 * i + 1])^6 +
    (tan(x[2 * i + 1] - x[2 * i + 2]) + x[2 * i + 1] - x[2 * i + 2])^4 +
    x[2 * i - 1]^8 +
    (x[2 * i + 2] - 1)^2 for i = 1:(div(n, 2) - 1)
  )
end
start_cragglvy(n::Int) = begin
  x0 = start_ones(n)
  x0[1] = 1
  return x0
end
cragglvy_ADNLPModel(n::Int = 100) =
  ADNLPModel(cragglvy, start_cragglvy(n), name = "cragglvy " * string(n) * " variables")

curly10(x::AbstractVector{Y}) where {Y <: Number} = curly(x, b = 10)
curly20(x::AbstractVector{Y}) where {Y <: Number} = curly(x, b = 20)
curly30(x::AbstractVector{Y}) where {Y <: Number} = curly(x, b = 30)
function curly(x::AbstractVector{Y}; b::Int = 10) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("curly: number of variables must be ≥ 2")
  n = max(2, n)

  f = Vector{Y}(undef, n)
  map!(i -> sum(x[j] for j = i:min(i + b, n)), f, [1:n;])

  sum(f[i] * (f[i] * (f[i]^2 - 20) - 0.1) for i = 1:n)
end
start_curly(n::Int) = [1.0e-4 * i / (n + 1) for i = 1:n]
curly_ADNLPModel(n::Int = 100) =
  ADNLPModel(curly, start_cragglvy(n), name = "curly " * string(n) * " variables")
curly10_ADNLPModel(n::Int = 100) =
  ADNLPModel(curly10, start_cragglvy(n), name = "curly10 " * string(n) * " variables")
curly20_ADNLPModel(n::Int = 100) =
  ADNLPModel(curly20, start_cragglvy(n), name = "curly20 " * string(n) * " variables")
curly30_ADNLPModel(n::Int = 100) =
  ADNLPModel(curly30, start_cragglvy(n), name = "curly30 " * string(n) * " variables")

dixmaanf(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaane(x, α = 1.0, β = 0.0625, γ = 0.0625, δ = 0.0625)
dixmaang(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaane(x, α = 1.0, β = 0.125, γ = 0.125, δ = 0.125)
dixmaanh(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaane(x, α = 1.0, β = 0.26, γ = 0.26, δ = 0.26)
function dixmaane(
  x::AbstractVector{Y};
  α::Float64 = 1.0,
  β::Float64 = 0.0,
  γ::Float64 = 0.125,
  δ::Float64 = 0.125,
) where {Y <: Number}
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum(i / n * α * x[i]^2 for i = 1:n) +
  sum(β * x[i]^2 * (x[i + 1] + x[i + 1]^2)^2 for i = 1:(n - 1)) +
  sum(γ * x[i]^2 * x[i + m]^4 for i = 1:(2 * m)) +
  sum(i / n * δ * x[i] * x[i + 2 * m] for i = 1:m)
end
start_dixmaane(n::Int) = [2.0 for i = 1:n]
dixmaane_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaane, start_dixmaane(n), name = "dixmaane " * string(n) * " variables")
dixmaanf_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaanf, start_dixmaane(n), name = "dixmaanf " * string(n) * " variables")
dixmaang_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaang, start_dixmaane(n), name = "dixmaang " * string(n) * " variables")
dixmaanh_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaanh, start_dixmaane(n), name = "dixmaanh " * string(n) * " variables")

dixmaanj(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaani(x, α = 1.0, β = 0.0625, γ = 0.0625, δ = 0.0625)
dixmaank(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaani(x, α = 1.0, β = 0.125, γ = 0.125, δ = 0.125)
dixmaanl(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaani(x, α = 1.0, β = 0.26, γ = 0.26, δ = 0.26)
function dixmaani(
  x::AbstractVector{Y};
  α::Float64 = 1.0,
  β::Float64 = 0.0,
  γ::Float64 = 0.125,
  δ::Float64 = 0.125,
) where {Y <: Number}
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum((i / n)^2 * α * x[i]^2 for i = 1:n) +
  sum(β * x[i]^2 * (x[i + 1] + x[i + 1]^2)^2 for i = 1:(n - 1)) +
  sum(γ * x[i]^2 * x[i + m]^4 for i = 1:(2 * m)) +
  sum((i / n)^2 * δ * x[i] * x[i + 2 * m] for i = 1:m)
end
start_dixmaani(n::Int) = [2.0 for i = 1:n]
dixmaani_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaani, start_dixmaani(n), name = "dixmaani " * string(n) * " variables")
dixmaanj_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaanj, start_dixmaani(n), name = "dixmaanj " * string(n) * " variables")
dixmaank_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaank, start_dixmaani(n), name = "dixmaank " * string(n) * " variables")
dixmaanl_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaanl, start_dixmaani(n), name = "dixmaanl " * string(n) * " variables")

dixmaann(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaanm(x, α = 1.0, β = 0.0625, γ = 0.0625, δ = 0.0625)
dixmaano(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaanm(x, α = 1.0, β = 0.125, γ = 0.125, δ = 0.125)
dixmaanp(x::AbstractVector{Y}) where {Y <: Number} =
  dixmaanm(x, α = 1.0, β = 0.26, γ = 0.26, δ = 0.26)
function dixmaanm(
  x::AbstractVector{Y};
  α::Float64 = 1.0,
  β::Float64 = 0.0,
  γ::Float64 = 0.125,
  δ::Float64 = 0.125,
) where {Y <: Number}
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum((i / n)^2 * α * x[i]^2 for i = 1:n) +
  sum(i / n * β * x[i]^2 * (x[i + 1] + x[i + 1]^2)^2 for i = 1:(n - 1)) +
  sum(i / n * γ * x[i]^2 * x[i + m]^4 for i = 1:(2 * m)) +
  sum((i / n)^2 * δ * x[i] * x[i + 2 * m] for i = 1:m)
end
start_dixmaanm(n::Int) = [2.0 for i = 1:n]
dixmaanm_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaanm, start_dixmaanm(n), name = "dixmaanm " * string(n) * " variables")
dixmaann_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaann, start_dixmaanm(n), name = "dixmaann " * string(n) * " variables")
dixmaano_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaano, start_dixmaanm(n), name = "dixmaano " * string(n) * " variables")
dixmaanp_ADNLPModel(n::Int = 99) =
  ADNLPModel(dixmaanp, start_dixmaanm(n), name = "dixmaanp " * string(n) * " variables")

function dixon3dq(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  (x[1] - 1.0)^2 + (x[n] - 1.0)^2 + sum((x[i] - x[i + 1])^2 for i = 2:(n - 1))
end
start_dixon3dq(n::Int) = start_ones(n)
dixon3dq_ADNLPModel(n::Int = 100) =
  ADNLPModel(dixon3dq, start_dixon3dq(n), name = "dixon3dq " * string(n) * " variables")

function dqdrtic(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)

  sum(x[i]^2 + 100 * (x[i + 1]^2 + x[i + 2]^2) for i = 1:(n - 2))
end
start_dqdrtic(n::Int) = (x -> 3 * x).(start_ones(n))
dqdrtic_ADNLPModel(n::Int = 100) =
  ADNLPModel(dqdrtic, start_dqdrtic(n), name = "dqdrtic " * string(n) * " variables")

function dqrtic(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)

  sum((x[i] - i)^4 for i = 1:n)
end
start_dqrtic(n::Int) = (x -> 2 * x).(start_ones(n))
dqrtic_ADNLPModel(n::Int = 100) =
  ADNLPModel(dqrtic, start_dqrtic(n), name = "dqrtic " * string(n) * " variables")

function edensch(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("edensch: number of variables must be ≥ 2")
  n = max(2, n)

  16 + sum((x[i] - 2)^4 + (x[i] * x[i + 1] - 2 * x[i + 1])^2 + (x[i + 1] + 1)^2 for i = 1:(n - 1))
end
start_edensch(n::Int) = (x -> 0 * x).(start_ones(n))
edensch_ADNLPModel(n::Int = 100) =
  ADNLPModel(edensch, start_edensch(n), name = "edensch " * string(n) * " variables")

function eg2(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("eg2: number of variables must be ≥ 2")
  n = max(2, n)

  sum(sin(x[1] + x[i]^2 - 1) for i = 1:(n - 1)) + 0.5 * sin(x[n]^2)
end
start_eg2(n::Int) = (x -> 0 * x).(start_ones(n))
eg2_ADNLPModel(n::Int = 100) =
  ADNLPModel(eg2, start_eg2(n), name = "eg2 " * string(n) * " variables")

function engval1(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("engval1: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i]^2 + x[i + 1]^2)^2 - 4 * x[i] + 3 for i = 1:(n - 1))
end
start_engval1(n::Int) = (x -> 2 * x).(start_ones(n))
engval1_ADNLPModel(n::Int = 100) =
  ADNLPModel(engval1, start_engval1(n), name = "engval1 " * string(n) * " variables")

function errinros_mod(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("errinros_mod: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i - 1] - 16.0 * x[i]^2 * (1.5 + sin(i))^2)^2 for i = 2:n) + sum((1.0 - x[i])^2 for i = 2:n)
end
start_errinros_mod(n::Int) = (x -> -1 * x).(start_ones(n))
errinros_mod_ADNLPModel(n::Int = 100) =
  ADNLPModel(errinros_mod, start_errinros_mod(n), name = "errinros_mod " * string(n) * " variables")

function freuroth(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("freuroth: number of variables must be ≥ 2")
  n = max(2, n)
  ngs = n - 1

  sum(((5.0 - x[i + 1]) * x[i + 1]^2 + x[i] - 2 * x[i + 1] - 13.0)^2 for i = 1:ngs) +
  sum(((1.0 + x[i + 1]) * x[i + 1]^2 + x[i] - 14 * x[i + 1] - 29.0)^2 for i = 1:ngs)
end
start_freuroth(n::Int) = begin
  x0 = zeros(n)
  x0[1] = 0.5
  x0[2] = -2.0
  return x0
end
freuroth_ADNLPModel(n::Int = 100) =
  ADNLPModel(freuroth, start_freuroth(n), name = "freuroth " * string(n) * " variables")

function genhumps(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)

  ζ = 20.0
  sum((sin(ζ * x[i])^2 * sin(ζ * x[i + 1])^2 + 0.05 * (x[i]^2 + x[i + 1]^2)) for i = 1:(n - 1))
end
start_genhumps(n::Int) = begin
  x0 = (x -> -506.2 * x).(start_ones(n))
  x0[1] = -506.0
  return x0
end
genhumps_ADNLPModel(n::Int = 100) =
  ADNLPModel(genhumps, start_genhumps(n), name = "genhumps " * string(n) * " variables")

function liarwhd(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("liarwhd: number of variables must be ≥ 4")
  n = max(2, n)

  sum(4.0 * (x[i]^2 - x[1])^2 + (x[i] - 1)^2 for i = 1:n)
end
start_liarwhd(n::Int) = (x -> 4 * x).(start_ones(n))
liarwhd_ADNLPModel(n::Int = 100) =
  ADNLPModel(liarwhd, start_liarwhd(n), name = "liarwhd " * string(n) * " variables")

function morebv(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("morebv: number of variables must be ≥ 4")
  n = max(2, n)

  h = 1.0 / (n + 1)

  sum((2.0 * x[i] - x[i - 1] - x[i + 1] + (h^2 / 2.0) * (x[i] + i * h + 1)^3)^2 for i = 2:(n - 1)) +
  (2.0 * x[1] - x[2] + (h^2 / 2.0) * (x[1] + 1)^3)^2 +
  (2.0 * x[n] - x[n - 1] + (h^2 / 2.0) * (x[n] + n * h + 1)^3)^2
end
start_morebv(n::Int) = (x -> 0.5 * x).(start_ones(n))
morebv_ADNLPModel(n::Int = 100) =
  ADNLPModel(morebv, start_morebv(n), name = "morebv " * string(n) * " variables")

function noncvxu2(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("noncvxu2: number of variables must be ≥ 2")
  n = max(2, n)

  sum(
    (x[i] + x[mod(3 * i - 2, n) + 1] + x[mod(7 * i - 3, n) + 1])^2 +
    4.0 * cos(x[i] + x[mod(3 * i - 2, n) + 1] + x[mod(7 * i - 3, n) + 1]) for i = 1:n
  )
end
start_noncvxu2(n::Int) = (x -> (Float64)(x)).([1:n;])
noncvxu2_ADNLPModel(n::Int = 100) =
  ADNLPModel(noncvxu2, start_noncvxu2(n), name = "noncvxu2 " * string(n) * " variables")

function noncvxun(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("noncvxun: number of variables must be ≥ 2")
  n = max(2, n)

  sum(
    (x[i] + x[mod(2 * i - 1, n) + 1] + x[mod(3 * i - 1, n) + 1])^2 +
    4.0 * cos(x[i] + x[mod(2 * i - 1, n) + 1] + x[mod(3 * i - 1, n) + 1]) for i = 1:n
  )
end
start_noncvxun(n::Int) = (x -> (Float64)(x)).([1:n;])
noncvxun_ADNLPModel(n::Int = 100) =
  ADNLPModel(noncvxun, start_noncvxun(n), name = "noncvxun " * string(n) * " variables")

function nondquar(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("nondquar: number of variables must be ≥ 2")
  n = max(2, n)

  (x[1] - x[2])^2 + (x[n - 1] - x[n])^2 + sum((x[i] + x[i + 1] + x[n])^4 for i = 1:(n - 2))
end
start_nondquar(n::Int) = begin
  x0 = ones(n)
  x0[2 * collect(1:div(n, 2))] .= -1.0
  return x0
end
nondquar_ADNLPModel(n::Int = 100) =
  ADNLPModel(nondquar, start_nondquar(n), name = "nondquar " * string(n) * " variables")

function quartc(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)

  sum((x[i] - i)^4 for i = 1:n)
end
start_quartc(n::Int) = (x -> 2 * x).(start_ones(n))
quartc_ADNLPModel(n::Int = 100) =
  ADNLPModel(quartc, start_quartc(n), name = "quartc " * string(n) * " variables")

function sbrybnd(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("sbrybnd: number of variables must be ≥ 2")
  n = max(2, n)

  p = zeros(n)
  J = Array{Any}(undef, n)
  for i = 1:n
    p[i] = exp(6.0 * (i - 1) / (n - 1))
    J[i] = [max(1, i - 5):(i - 1); (i + 1):min(n, i + 1)]
  end

  sum(
    (
      (2.0 + 5.0 * p[i]^2 * x[i]^2) * p[i] * x[i] + 1.0 -
      sum(p[j] * x[j] * (1.0 + p[j] * x[j]) for j in J[i])
    )^2 for i = 1:n
  )
end
function start_sbrybnd(n::Int)
  p = zeros(n)
  for i = 1:n
    p[i] = exp(6.0 * (i - 1) / (n - 1))
  end
  x0 = map(pᵢ -> 1.0 / pᵢ, p)
  return x0
end
sbrybnd_ADNLPModel(n::Int = 100) =
  ADNLPModel(sbrybnd, start_sbrybnd(n), name = "sbrybnd " * string(n) * " variables")

function tridia(
  x::AbstractVector{Y},
  α::Float64 = 2.0,
  β::Float64 = 1.0,
  γ::Float64 = 1.0,
  δ::Float64 = 1.0,
) where {Y <: Number}
  n = length(x)

  γ * (x[1] * δ - 1.0)^2 + sum(i * (-β * x[i - 1] + α * x[i])^2 for i = 2:n)
end
start_tridia(n::Int) = start_ones(n)
tridia_ADNLPModel(n::Int = 100) =
  ADNLPModel(tridia, start_tridia(n), name = "tridia " * string(n) * " variables")

function scosine(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("scosine: number of variables must be ≥ 2")
  n = max(2, n)

  p = zeros(n)
  for i = 1:n
    p[i] = exp(6.0 * (i - 1) / (n - 1))
  end

  sum(cos(p[i]^2 * x[i]^2 - p[i + 1] * x[i + 1] / 2.0) for i = 1:(n - 1))
end
function start_scosine(n::Int)
  p = zeros(n)
  for i = 1:n
    p[i] = exp(6.0 * (i - 1) / (n - 1))
  end
  x0 = map(pᵢ -> 1.0 / pᵢ, p)
  return x0
end
scosine_ADNLPModel(n::Int = 100) =
  ADNLPModel(scosine, start_scosine(n), name = "scosine " * string(n) * " variables")

function sinquad(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 3 && @warn("sinquad: number of variables must be ≥ 3")
  n = max(3, n)

  (x[1] - 1.0)^4 +
  (x[n]^2 - x[1]^2)^2 +
  sum((sin(x[i] - x[n]) - x[1]^2 + x[i]^2)^2 for i = 2:(n - 1))
end
start_sinquad(n::Int) = (x -> 0.1 * x).(start_ones(n))
sinquad_ADNLPModel(n::Int = 100) =
  ADNLPModel(sinquad, start_sinquad(n), name = "sinquad " * string(n) * " variables")

function srosenbr(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  (n % 2 == 0) || @warn("srosenbr: number of variables adjusted to be even")
  n = 2 * max(1, div(n, 2))

  sum(100.0 * (x[2 * i] - x[2 * i - 1]^2)^2 + (x[2 * i - 1] - 1.0)^2 for i = 1:div(n, 2))
end
start_srosenbr(n::Int) = begin
  x0 = ones(n)
  x0[2 * (collect(1:div(n, 2))) .- 1] .= -1.2
  return x0
end
srosenbr_ADNLPModel(n::Int = 100) =
  ADNLPModel(srosenbr, start_srosenbr(n), name = "srosenbr " * string(n) * " variables")

function woods(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  (n % 4 == 0) || @warn("woods: number of variables adjusted to be a multiple of 4")
  n = 4 * max(1, div(n, 4))

  1.0 + sum(
    100 * (x[4 * i - 2] - x[4 * i - 3]^2)^2 +
    (1 - x[4 * i - 3])^2 +
    90 * (x[4 * i] - x[4 * i - 1]^2)^2 +
    (1 - x[4 * i - 1])^2 +
    10 * (x[4 * i - 2] + x[4 * i] - 2)^2 +
    0.1 * (x[4 * i - 2] - x[4 * i])^2 for i = 1:div(n, 4)
  )
end

start_woods(n::Int) = (x -> -2 * x).(start_ones(n))
woods_ADNLPModel(n::Int = 100) =
  ADNLPModel(woods, start_woods(n), name = "woods " * string(n) * " variables")

# test(x :: AbstractVector{Y}) where Y <: Number =n = length(x)
# start_(n :: Int) =  (x -> 2 * x).(start_ones(n))
# _ADNLPModel(n :: Int=100) = ADNLPModel(eg2, start_(n))

# name_arwhead(x :: AbstractVector{Y}) where Y <: Number = "arwhead" * string(length(x))
# arwhead_ADNLPModel(x :: AbstractVector{Y}) where Y <: Number = ADNLPModel(eval(arwhead), name=name_arwhead(x))

#
# fx = arwhead(x)
# @show fx
# x = ones(50)
# n = 52
# dixmaanf_adnlp = woods_ADNLPModel(n)
# ges = JSOSolvers.trunk(dixmaanf_adnlp)
# @show ges

#peu proprice par manque de prétraitement de l'arbre
function extrosnb(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("extrosnb: number of variables must be ≥ 2")
  n = max(2, n)

  100.0 * sum((x[i] - x[i - 1]^2)^2 for i = 2:n) + (1.0 - x[1])^2
end
start_extrosnb(n::Int) = (x -> -1 * x).(start_ones(n))
extrosnb_ADNLPModel(n::Int = 100) = ADNLPModel(extrosnb, start_extrosnb(n))

function chnrosnb_mod(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("chnrosnb: number of variables must be ≥ 2")
  n = max(2, n)

  16 * sum((x[i - 1] - x[i]^2)^2 * (1.5 + sin(i))^2 for i = 2:n) + sum((1.0 - x[i])^2 for i = 2:n)
end
start_chnrosnb_mod(n::Int) = start_ones(n)
chnrosnb_mod_ADNLPModel(n::Int = 100) = ADNLPModel(chnrosnb_mod, start_chnrosnb_mod(n))

function nondia(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 2 && @warn("nondia: number of variables must be ≥ 2")
  n = max(2, n)

  (x[1] - 1.0)^2 + 100 * sum((x[1] - x[i]^2)^2 for i = 2:n)
end
start_nondia(n::Int) = (x -> -1 * x).(start_ones(n))
nondia_ADNLPModel(n::Int = 100) = ADNLPModel(nondia, start_nondia(n))

# pas de SPS
function power(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)

  (sum((i * x[i]^2) for i = 1:n))^2
end
start_power(n::Int) = start_ones(n)
power_ADNLPModel(n::Int = 100) = ADNLPModel(power, start_power(n))

function vardim(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  nlp = Model()

  sum((x[i] - 1)^2 for i = 1:n) +
  (sum(i * (x[i] - 1) for i = 1:n))^2 +
  (sum(i * (x[i] - 1) for i = 1:n))^4
end
start_vardim(n::Int) = map(i -> (Float64)(1 - i / n), [1:n;])
vardim_ADNLPModel(n::Int = 100) = ADNLPModel(vardim, start_vardim(n))

function sparsine(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 10 && @warn("sparsine: number of variables must be ≥ 10")
  n = max(10, n)

  0.5 * sum(
    i *
    (
      sin(x[i]) +
      sin(x[mod(2 * i - 1, n) + 1]) +
      sin(x[mod(3 * i - 1, n) + 1]) +
      sin(x[mod(5 * i - 1, n) + 1]) +
      sin(x[mod(7 * i - 1, n) + 1]) +
      sin(x[mod(11 * i - 1, n) + 1])
    )^2 for i = 1:n
  )
end
start_sparsine(n::Int) = (x -> 0.5 * x).(start_ones(n))
sparsine_ADNLPModel(n::Int = 100) = ADNLPModel(sparsine, start_sparsine(n))

function sparsqur(x::AbstractVector{Y}) where {Y <: Number}
  n = length(x)
  n < 10 && @warn("sparsqur: number of variables must be ≥ 10")
  n = max(10, n)

  1 / 8 * sum(
    i *
    (
      x[i]^2 +
      x[mod(2 * i - 1, n) + 1]^2 +
      x[mod(3 * i - 1, n) + 1]^2 +
      x[mod(5 * i - 1, n) + 1]^2 +
      x[mod(7 * i - 1, n) + 1]^2 +
      x[mod(11 * i - 1, n) + 1]^2
    )^2 for i = 1:n
  )
end
start_sparsqur(n::Int) = (x -> 0.5 * x).(start_ones(n))
sparsqur_ADNLPModel(n::Int = 100) = ADNLPModel(sparsqur, start_sparsqur(n))
