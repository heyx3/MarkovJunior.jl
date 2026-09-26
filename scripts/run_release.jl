using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
println(stderr, ".")

using MarkovJunior
markovjunior_run_gui()