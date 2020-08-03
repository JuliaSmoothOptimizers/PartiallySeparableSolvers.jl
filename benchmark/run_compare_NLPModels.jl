using PkgBenchmark
using SolverBenchmark
using Plots

ENV["GKSwstype"]=100

# perform benchmarks
results = PkgBenchmark.benchmarkpkg("PartiallySeparableSolvers", script="benchmark/scripts/compare_with_NLPModels.jl")
p = profile_solvers(results)
savefig(p, "benchmark/profiles/profile_NLPModels.pdf")


results_obj = PkgBenchmark.benchmarkpkg("PartiallySeparableSolvers", script="benchmark/scripts/compare_obj.jl")
p_obj = profile_solvers(results_obj)
savefig(p_obj, "benchmark/profiles/profile_obj.pdf")

results_grad = PkgBenchmark.benchmarkpkg("PartiallySeparableSolvers", script="benchmark/scripts/compare_grad.jl")
p_grad = profile_solvers(results_grad)
savefig(p_grad, "benchmark/profiles/profile_grad.pdf")

results_Hv = PkgBenchmark.benchmarkpkg("PartiallySeparableSolvers", script="benchmark/scripts/compare_Hv.jl")
p_Hv = profile_solvers(results_Hv)
savefig(p_Hv, "benchmark/profiles/profile_Hv.pdf")
