using PkgBenchmark
using SolverBenchmark
using PartiallySeparableSolvers
using Plots


commit = benchmarkpkg("PartiallySeparableSolvers",script="benchmark/solver_benchmark.jl")  #dernier commit sur la branche sur laquelle on se trouve
master = benchmarkpkg("PartiallySeparableSolvers", "master", script="benchmark/solver_benchmark.jl") # branche masterjudgement_solver = judge(master_solver, commit_solver)
judgement_solver = judge(master, commit)

export_markdown("benchmark/judgement_solver.md", judgement_solver)
ENV["GKSwstype"]=100
p = SolverBenchmark.profile_solvers(commit)
savefig(p, "benchmark/profile_solver.pdf")
