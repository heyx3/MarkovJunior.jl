# Pass "-loglogic" to turn on algorithm logging.
# Pass "-logmd" to turn on multidimensional-rewrite rule logging.

using Pkg
Pkg.activate(joinpath(@__FILE__, "../.."))
println(".")

using MarkovJunior
MarkovJunior.markovjunior_asserts_enabled() = true
if "-loglogic" in ARGS
    MarkovJunior.log_logic() = true
end
if "-logmd" in ARGS
    MarkovJunior.log_md_symmetry_logic() = true
end
markovjunior_run_tool()