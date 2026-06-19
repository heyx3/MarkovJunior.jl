# Runs the IPC process without having to compile this project into an executable,
#    for development purposes.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
println(stderr, ".")

using MarkovJunior
MarkovJunior.markovjunior_run_ipc_main()