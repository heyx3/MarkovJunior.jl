# Make sure the test is always running in the same directory and within the same project.
using Pkg
Pkg.activate(@__DIR__)
insert!(LOAD_PATH, 1, joinpath(@__DIR__, ".."))

using Sockets, UUIDs
using Suppressor # Capture IPC server's stderr

using MarkovJunior; const MJ = MarkovJunior
MJ.markovjunior_asserts_enabled() = true

using Bplus; @using_bplus
@bp_check(Bplus.BplusCore === MJ.Bplus.BplusCore,
          "Test project isn't using the same version of B+")

###################################################


include("utils.jl")
include("parsing.jl")
include("md_symmetries.jl")
include("add_ons.jl")
include("ipc.jl")

println("\n\nTests passed!\n")