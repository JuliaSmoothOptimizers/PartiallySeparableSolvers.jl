using JuMP, ADNLPModels, JSOSolvers



start_ones(n :: Int) = ones(n)

function arwhead(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("arwhead: number of variables must be ≥ 2")
  n = max(2, n)

  return sum((x[i]^2 + x[n]^2)^2 - 4 * x[i] + 3 for i=1:n-1)
end
start_arwhead(n :: Int) = ones(n)
arwhead_ADNLPModel(n :: Int=100) = RADNLPModel(arwhead, start_arwhead(n))



function bdqrtic(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 5 && @warn("bdqrtic: number of variables must be ≥ 5")
  n = max(5, n)

  return sum((3 - 4 * x[i])^2 + (x[i]^2 + 2 * x[i+1]^2 + 3 * x[i+2]^2 + 4 * x[i+3]^2 + 5 * x[n]^2)^2 for i=1:n-4)
end
start_bdqrtic(n :: Int) = ones(n)
bdqrtic_ADNLPModel(n :: Int=100) where Y <: Number = RADNLPModel(bdqrtic, start_bdqrtic(n))


function brybnd(x :: AbstractVector{Y}; ml :: Int=5, mu :: Int=1) where Y <: Number
  n = length(x)

  sum(
    (
      x[i] * (2 + 5 * x[i]^2) + 1 -
      sum(
        x[j] * (1 + x[j])
        for j = max(1, i-ml) : min(n, i+mu) if j != i
      )
    )^2 for i=1:n
  )
end
start_brybnd(n :: Int) = (x -> -1 * x).(ones(n))
brybnd_ADNLPModel(n :: Int=100) = RADNLPModel(brybnd, start_bdqrtic(n))


function chainwoo(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  (n % 4 == 0) || @warn("chainwoo: number of variables adjusted to be a multiple of 4")
  n = 4 * max(1, div(n, 4))

  1.0 + sum(100 * (x[2*i]   - x[2*i-1]^2)^2 + (1 - x[2*i-1])^2 +
               90 * (x[2*i+2] - x[2*i+1]^2)^2 + (1 - x[2*i+1])^2 +
               10 * (x[2*i] + x[2*i+2] - 2)^2 + 0.1 * (x[2*i] - x[2*i+2])^2 for i=1:div(n,2)-1)
end
function start_chainwoo(n :: Int)
  x0 = (x -> -2 * x).(ones(n))
  x0[1] = -3
  x0[2] = -1
  x0[3] = -3
  x0[4] = -1
  return x0
end
chainwoo_ADNLPModel(n :: Int=100) = RADNLPModel(chainwoo, start_chainwoo(n))


function chnrosnb_mod(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("chnrosnb: number of variables must be ≥ 2")
  n = max(2, n)

  16 * sum((x[i-1] - x[i]^2)^2*(1.5+sin(i))^2 for i=2:n) + sum((1.0 - x[i])^2 for i=2:n)
end
start_chnrosnb_mod(n :: Int) = start_ones(n)
chnrosnb_mod_ADNLPModel(n :: Int=100) = RADNLPModel(chnrosnb_mod, start_chnrosnb_mod(n))


function cosine(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("cosine: number of variables must be ≥ 2")
  n = max(2, n)
  sum(cos(x[i]^2 - 0.5 * x[i+1]) for i = 1:n-1)
end
start_cosine(n :: Int) = start_ones(n)
cosine_ADNLPModel(n :: Int=100) = RADNLPModel(cosine, start_cosine(n))


function cragglvy(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("cragglvy: number of variables must be ≥ 2")
  n = max(2, n)

  sum((exp(x[2*i-1]) - x[2*i])^4 + 100 * (x[2*i] - x[2*i+1])^6 +
        (tan(x[2*i+1] - x[2*i+2]) + x[2*i+1] - x[2*i+2])^4 +
        x[2*i-1]^8 + (x[2*i+2] - 1)^2 for i = 1:div(n,2)-1)
end
start_cragglvy(n :: Int) = begin x0 = start_ones(n); x0[1] = 1; return x0 end
cragglvy_ADNLPModel(n :: Int=100) = RADNLPModel(cragglvy, start_cragglvy(n))


curly10(x :: AbstractVector{Y}) where Y <: Number = curly(x, b=10)
curly20(x :: AbstractVector{Y}) where Y <: Number = curly(x, b=20)
curly30(x :: AbstractVector{Y}) where Y <: Number = curly(x, b=30)
function curly(x :: AbstractVector{Y}; b :: Int=10) where Y <: Number
  n = length(x)
  n < 2 && @warn("curly: number of variables must be ≥ 2")
  n = max(2, n)

  f = Vector{Y}(undef,n)
  map!(i -> sum(x[j] for j=i:min(i+b,n)),f, [1:n;])

  sum(f[i] * (f[i] * (f[i]^2 - 20) - 0.1) for i = 1:n)
end
start_curly(n :: Int) = [1.0e-4 * i /(n+1) for i = 1:n]
curly_ADNLPModel(n :: Int=100) = RADNLPModel(curly, start_cragglvy(n))
curly10_ADNLPModel(n :: Int=100) = RADNLPModel(curly10, start_cragglvy(n))
curly20_ADNLPModel(n :: Int=100) = RADNLPModel(curly20, start_cragglvy(n))
curly30_ADNLPModel(n :: Int=100) = RADNLPModel(curly30, start_cragglvy(n))


dixmaanf(x :: AbstractVector{Y}) where Y <: Number = dixmaane(x, α=1.0, β=0.0625, γ=0.0625, δ=0.0625)
dixmaang(x :: AbstractVector{Y}) where Y <: Number = dixmaane(x, α=1.0, β=0.125, γ=0.125, δ=0.125)
dixmaanh(x :: AbstractVector{Y}) where Y <: Number = dixmaane(x, α=1.0, β=0.26, γ=0.26, δ=0.26)
function dixmaane(x :: AbstractVector{Y};
                  α :: Float64=1.0, β :: Float64=0.0, γ :: Float64=0.125, δ :: Float64=0.125) where Y <: Number
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum(i / n * α * x[i]^2 for                 i=1:n) +
  sum(β * x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i=1:n-1) +
  sum(γ * x[i]^2 * x[i+m]^4 for              i=1:2*m) +
  sum(i / n * δ * x[i] * x[i+2*m] for        i=1:m)
end
start_dixmaane(n :: Int) = [2.0 for i = 1:n]
dixmaane_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaane, start_dixmaane(n))
dixmaanf_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanf, start_dixmaane(n))
dixmaang_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaang, start_dixmaane(n))
dixmaanh_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanh, start_dixmaane(n))


dixmaanj(x :: AbstractVector{Y}) where Y <: Number = dixmaani(x, α=1.0, β=0.0625, γ=0.0625, δ=0.0625)
dixmaank(x :: AbstractVector{Y}) where Y <: Number = dixmaani(x, α=1.0, β=0.125, γ=0.125, δ=0.125)
dixmaanl(x :: AbstractVector{Y}) where Y <: Number = dixmaani(x, α=1.0, β=0.26, γ=0.26, δ=0.26)
function dixmaani(x :: AbstractVector{Y};
                  α :: Float64=1.0, β :: Float64=0.0, γ :: Float64=0.125, δ :: Float64=0.125) where Y <: Number
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum((i / n)^2 * α * x[i]^2 for             i=1:n)   +
  sum(β * x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i=1:n-1) +
  sum(γ * x[i]^2 * x[i+m]^4 for              i=1:2*m) +
  sum((i / n)^2 * δ * x[i] * x[i+2*m] for    i=1:m)
end
start_dixmaani(n :: Int) = [2.0 for i = 1:n]
dixmaani_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaani, start_dixmaani(n))
dixmaanj_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanj, start_dixmaani(n))
dixmaank_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaank, start_dixmaani(n))
dixmaanl_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanl, start_dixmaani(n))


dixmaann(x :: AbstractVector{Y}) where Y <: Number = dixmaanm(x, α=1.0, β=0.0625, γ=0.0625, δ=0.0625)
dixmaano(x :: AbstractVector{Y}) where Y <: Number = dixmaanm(x, α=1.0, β=0.125, γ=0.125, δ=0.125)
dixmaanp(x :: AbstractVector{Y}) where Y <: Number = dixmaanm(x, α=1.0, β=0.26, γ=0.26, δ=0.26)
function dixmaanm(x :: AbstractVector{Y};
                  α :: Float64=1.0, β :: Float64=0.0, γ :: Float64=0.125, δ :: Float64=0.125) where Y <: Number
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum((i / n)^2 * α * x[i]^2 for                     i=1:n) +
  sum(i / n * β * x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i=1:n-1) +
  sum(i / n * γ * x[i]^2 * x[i+m]^4 for              i=1:2*m) +
  sum((i / n)^2 * δ * x[i] * x[i+2*m] for            i=1:m)

end
start_dixmaanm(n :: Int) = [2.0 for i = 1:n]
dixmaanm_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanm, start_dixmaanm(n))
dixmaann_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaann, start_dixmaanm(n))
dixmaano_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaano, start_dixmaanm(n))
dixmaanp_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanp, start_dixmaanm(n))



function dixon3dq(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  (x[1] - 1.0)^2 + (x[n] - 1.0)^2 + sum((x[i] - x[i+1])^2 for i=2:n-1)
end
start_dixon3dq(n :: Int) =  start_ones(n)
dixon3dq_ADNLPModel(n :: Int=100) = RADNLPModel(dixon3dq, start_dixon3dq(n))



function dqdrtic(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  sum(x[i]^2 + 100 * (x[i+1]^2 + x[i+2]^2) for i=1:n-2)
end
start_dqdrtic(n :: Int) =  (x -> 3 * x).(start_ones(n))
dqdrtic_ADNLPModel(n :: Int=100) = RADNLPModel(dqdrtic, start_dqdrtic(n))


function dqrtic(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  sum((x[i] - i)^4 for i=1:n)
end
start_dqrtic(n :: Int) =  (x -> 2 * x).(start_ones(n))
dqrtic_ADNLPModel(n :: Int=100) = RADNLPModel(dqrtic, start_dqrtic(n))


function edensch(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("edensch: number of variables must be ≥ 2")
  n = max(2, n)

  16 +
  sum(
    (x[i] - 2)^4 + (x[i] * x[i+1] - 2 * x[i+1])^2 + (x[i+1] + 1)^2
    for i=1:n-1
  )
end
start_edensch(n :: Int) =  (x -> 0 * x).(start_ones(n))
edensch_ADNLPModel(n :: Int=100) = RADNLPModel(edensch, start_edensch(n))


function eg2(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("eg2: number of variables must be ≥ 2")
  n = max(2, n)

    sum(
      sin(x[1] + x[i]^2 - 1)
      for i=1:n-1
    ) +
    0.5 * sin(x[n]^2)
end
start_eg2(n :: Int) =  (x -> 0 * x).(start_ones(n))
eg2_ADNLPModel(n :: Int=100) = RADNLPModel(eg2, start_eg2(n))

function engval1(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("engval1: number of variables must be ≥ 2")
  n = max(2, n)

    sum(
      (x[i]^2 + x[i+1]^2)^2 - 4 * x[i] + 3
      for i=1:n-1
    )
end
start_engval1(n :: Int) =  (x -> 2 * x).(start_ones(n))
engval1_ADNLPModel(n :: Int=100) = RADNLPModel(engval1, start_engval1(n))



function errinros_mod(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("errinros_mod: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i-1] - 16.0 * x[i]^2 * (1.5 + sin(i))^2)^2 for i=2:n) + sum((1.0 - x[i])^2 for i=2:n)
end
start_errinros_mod(n :: Int) =  (x -> -1 * x).(start_ones(n))
errinros_mod_ADNLPModel(n :: Int=100) = RADNLPModel(errinros_mod, start_errinros_mod(n))


function extrosnb(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("extrosnb: number of variables must be ≥ 2")
  n = max(2, n)

  100.0 * sum((x[i] - x[i - 1]^2)^2 for i=2:n) + (1.0 - x[1])^2
end
start_extrosnb(n :: Int) =  (x -> -1 * x).(start_ones(n))
extrosnb_ADNLPModel(n :: Int=100) = RADNLPModel(extrosnb, start_extrosnb(n))


function freuroth(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("freuroth: number of variables must be ≥ 2")
  n = max(2, n)
  ngs = n - 1


  sum(((5.0 - x[i+1]) * x[i+1]^2 + x[i] - 2 * x[i+1] - 13.0)^2 for i=1:ngs) +
  sum(((1.0 + x[i+1]) * x[i+1]^2 + x[i] - 14 * x[i+1] - 29.0)^2 for i=1:ngs)
end
start_freuroth(n :: Int) = begin x0 = zeros(n); x0[1] = 0.5; x0[2] = -2.0; return x0 end
freuroth_ADNLPModel(n :: Int=100) = RADNLPModel(freuroth, start_freuroth(n))


function genhumps(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  ζ = 20.0
  sum(( sin(ζ * x[i])^2 * sin(ζ * x[i+1])^2 + 0.05 * (x[i]^2 + x[i+1]^2)) for i=1:n-1)
end
start_genhumps(n :: Int) =  begin x0 = (x -> -506.2 * x).(start_ones(n)); x0[1] = - 506.0; return x0 end
genhumps_ADNLPModel(n :: Int=100) = RADNLPModel(genhumps, start_genhumps(n))


function liarwhd(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("liarwhd: number of variables must be ≥ 4")
  n = max(2, n)

  sum(4.0*(x[i]^2 - x[1])^2 + (x[i] - 1)^2  for i=1:n)
end
start_liarwhd(n :: Int) =  (x -> 4 * x).(start_ones(n))
liarwhd_ADNLPModel(n :: Int=100) = RADNLPModel(liarwhd, start_liarwhd(n))


function morebv(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("morebv: number of variables must be ≥ 4")
  n = max(2, n)

  h = 1.0/(n+1)

  sum((2.0 * x[i] - x[i-1] - x[i+1] + (h^2 / 2.0) * (x[i] + i * h + 1)^3)^2 for i=2:n-1) +
  (2.0 * x[1] - x[2] + (h^2 / 2.0) * (x[1] + 1)^3)^2 +
  (2.0 * x[n] - x[n-1] + (h^2 / 2.0) * (x[n] + n * h + 1)^3)^2

end
start_morebv(n :: Int) =  (x -> 0.5 * x).(start_ones(n))
morebv_ADNLPModel(n :: Int=100) = RADNLPModel(morebv, start_morebv(n))

function noncvxu2(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("noncvxu2: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i] + x[mod(3 * i - 2, n) + 1] + x[mod(7 * i - 3, n) + 1])^2 +
  4.0 * cos(x[i] + x[mod(3 * i - 2, n) + 1] + x[mod(7 * i - 3, n) + 1]) for i=1:n)
end
start_noncvxu2(n :: Int) =  (x -> (Float64)(x)).([1:n;])
noncvxu2_ADNLPModel(n :: Int=100) = RADNLPModel(noncvxu2, start_noncvxu2(n))


function noncvxun(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("noncvxun: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i] + x[mod(2*i-1, n) + 1] + x[mod(3*i-1, n) + 1])^2 +
  4.0 * cos(x[i] + x[mod(2*i-1, n) + 1] + x[mod(3*i-1, n) + 1]) for i=1:n)
end
start_noncvxun(n :: Int) =  (x -> (Float64)(x)).([1:n;])
noncvxun_ADNLPModel(n :: Int=100) = RADNLPModel(noncvxun, start_noncvxun(n))


function nondia(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("nondia: number of variables must be ≥ 2")
  n = max(2, n)

  (x[1] - 1.0)^2 + 100*sum((x[1] - x[i]^2)^2 for i=2:n)
end
start_nondia(n :: Int) = (x -> -1 * x).(start_ones(n))
nondia_ADNLPModel(n :: Int=100) = RADNLPModel(nondia, start_nondia(n))


function nondquar(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("nondquar: number of variables must be ≥ 2")
  n = max(2, n)

  (x[1] - x[2])^2 + (x[n-1] - x[n])^2 + sum((x[i] + x[i+1] + x[n])^4 for i=1:n-2)
end
start_nondquar(n :: Int) = begin x0 = ones(n); x0[2 * collect(1:div(n, 2))] .= -1.0; return x0 end
nondquar_ADNLPModel(n :: Int=100) = RADNLPModel(nondquar, start_nondquar(n))

#pas partiellement séparable
function power(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  (sum((i * x[i]^2) for i=1:n))^2
end
start_power(n :: Int) =  start_ones(n)
power_ADNLPModel(n :: Int=100) = RADNLPModel(power, start_power(n))


function quartc(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  sum((x[i] - i)^4 for i=1:n)
end
start_quartc(n :: Int) =  (x -> 2 * x).(start_ones(n))
quartc_ADNLPModel(n :: Int=100) = RADNLPModel(quartc, start_quartc(n))


function sbrybnd(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("sbrybnd: number of variables must be ≥ 2")
  n = max(2, n)

  p = zeros(n)
  J = Array{Any}(undef, n)
  for i=1:n
    p[i] = exp(6.0*(i-1)/(n-1))
    J[i] = [max(1, i-5):i-1; i+1:min(n, i+1)]
  end

  sum(((2.0 + 5.0 * p[i]^2 * x[i]^2) * p[i] * x[i] + 1.0 - sum(p[j] * x[j] * (1.0 + p[j] * x[j]) for j=J[i]))^2 for i=1:n)
end
function start_sbrybnd(n :: Int)
  p = zeros(n)
  for i=1:n
    p[i] = exp(6.0*(i-1)/(n-1))
  end
  x0 = map(pᵢ -> 1.0/pᵢ, p)
  return x0
end
sbrybnd_ADNLPModel(n :: Int=100) = RADNLPModel(sbrybnd, start_sbrybnd(n))

function tridia(x :: AbstractVector{Y}, α::Float64=2.0, β::Float64=1.0, γ::Float64=1.0, δ::Float64=1.0) where Y <: Number
  n = length(x)

  γ * (x[1] * δ - 1.0)^2 + sum(i * (-β * x[i-1] + α * x[i])^2 for i=2:n)
end
start_tridia(n :: Int) = start_ones(n)
tridia_ADNLPModel(n :: Int=100) = RADNLPModel(tridia, start_tridia(n))


function vardim(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  nlp = Model()

  sum((x[i] - 1)^2 for i=1:n) + (sum(i * (x[i] - 1) for i=1:n))^2 + (sum(i * (x[i] - 1) for i=1:n))^4
end
start_vardim(n :: Int) =  map(i -> (Float64)(1 - i/n), [1:n;])
vardim_ADNLPModel(n :: Int=100) = RADNLPModel(vardim, start_vardim(n))


function scosine(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("scosine: number of variables must be ≥ 2")
  n = max(2, n)

  p = zeros(n)
  for i=1:n
    p[i] = exp(6.0 * (i-1) / (n-1))
  end

  sum(cos(p[i]^2 * x[i]^2 - p[i+1] * x[i+1] / 2.0) for i=1:n-1)
end
function start_scosine(n :: Int)
  p = zeros(n)
  for i=1:n
    p[i] = exp(6.0 * (i-1) / (n-1))
  end
  x0 = map(pᵢ -> 1.0/pᵢ, p)
  return x0
end
scosine_ADNLPModel(n :: Int=100) = RADNLPModel(scosine, start_scosine(n))


function sinquad(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 3 && @warn("sinquad: number of variables must be ≥ 3")
  n = max(3, n)

  (x[1] - 1.0)^4 + (x[n]^2 - x[1]^2)^2 + sum((sin(x[i] - x[n]) - x[1]^2 + x[i]^2)^2 for i=2:n-1)
end
start_sinquad(n :: Int) =  (x -> 0.1 * x).(start_ones(n))
sinquad_ADNLPModel(n :: Int=100) = RADNLPModel(sinquad, start_sinquad(n))


function sparsine(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 10 && @warn("sparsine: number of variables must be ≥ 10")
  n = max(10, n)

  0.5 * sum(
    i * (sin(x[i]) +
    sin(x[mod(2*i-1, n) + 1]) +
    sin(x[mod(3*i-1, n) + 1]) +
    sin(x[mod(5*i-1, n) + 1]) +
    sin(x[mod(7*i-1, n) + 1]) +
    sin(x[mod(11*i-1, n) + 1]))^2 for i=1:n)
end
start_sparsine(n :: Int) =  (x -> 0.5 * x).(start_ones(n))
sparsine_ADNLPModel(n :: Int=100) = RADNLPModel(sparsine, start_sparsine(n))


function sparsqur(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 10 && @warn("sparsqur: number of variables must be ≥ 10")
  n = max(10, n)

  1/8 * sum(
  i * (x[i]^2 +
  x[mod(2*i-1, n) + 1]^2 +
  x[mod(3*i-1, n) + 1]^2 +
  x[mod(5*i-1, n) + 1]^2 +
  x[mod(7*i-1, n) + 1]^2 +
  x[mod(11*i-1, n) + 1]^2)^2 for i=1:n)
end
start_sparsqur(n :: Int) =  (x -> 0.5 * x).(start_ones(n))
sparsqur_ADNLPModel(n :: Int=100) = RADNLPModel(sparsqur, start_sparsqur(n))


function srosenbr(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  (n % 2 == 0) || @warn("srosenbr: number of variables adjusted to be even")
  n = 2 * max(1, div(n, 2))

  sum(100.0 * (x[2*i] - x[2*i-1]^2)^2  + (x[2*i-1] - 1.0)^2 for i=1:div(n, 2))
end
start_srosenbr(n :: Int) = begin x0 = ones(n); x0[2*(collect(1:div(n,2))).-1] .= -1.2; return x0 end
srosenbr_ADNLPModel(n :: Int=100) = RADNLPModel(srosenbr, start_srosenbr(n))


function woods(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  (n % 4 == 0) || @warn("woods: number of variables adjusted to be a multiple of 4")
  n = 4 * max(1, div(n, 4))

  1.0 + sum(
    100 * (x[4*i-2] - x[4*i-3]^2)^2 + (1 - x[4*i-3])^2 +
    90 * (x[4*i] - x[4*i-1]^2)^2 + (1 - x[4*i-1])^2 +
    10 * (x[4*i-2] + x[4*i] - 2)^2 + 0.1 * (x[4*i-2] - x[4*i])^2 for i=1:div(n,4))
end

start_woods(n :: Int) =  (x -> -2 * x).(start_ones(n))
woods_ADNLPModel(n :: Int=100) = RADNLPModel(woods, start_woods(n))







test(x :: AbstractVector{Y}) where Y <: Number =n = length(x)
start_(n :: Int) =  (x -> 2 * x).(start_ones(n))
_ADNLPModel(n :: Int=100) = RADNLPModel(eg2, start_(n))

# name_arwhead(x :: AbstractVector{Y}) where Y <: Number = "arwhead" * string(length(x))
# arwhead_ADNLPModel(x :: AbstractVector{Y}) where Y <: Number = RADNLPModel(eval(arwhead), name=name_arwhead(x))

#
# fx = arwhead(x)
# @show fx
# x = ones(50)
# n = 52
# dixmaanf_adnlp = woods_ADNLPModel(n)
# ges = JSOSolvers.trunk(dixmaanf_adnlp)
# @show ges
