using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))
println(".")

using MarkovJunior
markovjunior_run_tool()