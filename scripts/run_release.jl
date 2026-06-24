using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))
println(stderr, ".")

using MarkovJunior
markovjunior_run_gui()