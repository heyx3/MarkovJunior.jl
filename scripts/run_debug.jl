using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))

using MarkovJunior
MarkovJunior.markovjunior_asserts_enabled() = true
MarkovJunior.main()