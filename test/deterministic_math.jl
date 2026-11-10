# Tests for `DeterministicMath`: accuracy, special cases, and identical bits across CPU targets.

const DM = MJ.DeterministicMath

(() -> begin

# The self-check passes (it already ran while precompiling and loading, or we wouldn't get here).
DM.check_determinism()

# Accuracy against Base, which is itself within a couple of ulps of the truth.
ulps(a, b) = (a == b) ? 0.0 : abs(a - b) / eps(abs(b))
spread(n, lo, hi) = [ lo + ((hi - lo) * (i / n)) for i in 1:n ]
function check_close(det_f, base_f, inputs, max_ulps)
    for x in inputs
        u = ulps(det_f(x...), base_f(x...))
        @bp_check(u <= max_ulps, det_f, x, " is ", u, " ulps from Base")
    end
end
check_close(DM.det_exp, exp, spread(10_000, -700, 700), 2)
check_close(DM.det_exp2, exp2, spread(10_000, -1000, 1000), 2)
check_close(DM.det_log, log, spread(10_000, 1e-300, 1e300), 2)
check_close(DM.det_log, log, spread(10_000, 0.5, 2.0), 4)
check_close(DM.det_log2, log2, spread(10_000, 1e-6, 1e6), 2)
check_close(DM.det_log10, log10, spread(10_000, 1e-6, 1e6), 2)
check_close(DM.det_pow, ^, [ (a, b) for a in spread(100, 0.001, 10.0) for b in spread(100, -8.0, 8.0) ], 40)
check_close(DM.det_pow, ^, [ (a, Float64(n)) for a in spread(100, -3.0, 3.0) for n in -20:20 ], 40)
check_close(DM.det_pow, ^, [ (Float32(a), Float32(b)) for a in spread(100, 0.001, 1.0) for b in spread(100, 0.3, 1.0) ], 1)

# Special cases behave like Base.
outcome(f, args...) = try f(args...) catch e; typeof(e) end
for (a, b) in ((0.0, 2.0), (0.0, -1.0), (-0.0, 3.0), (Inf, 2.0), (Inf, -2.0), (0.5, Inf), (2.0, -Inf),
               (-2.0, 3.0), (-2.0, 2.5), (NaN, 0.0), (1.0, NaN), (-8.0, 1e20), (2.0, 1100.0), (2.0, -1100.0))
    @bp_check(isequal(outcome(DM.det_pow, a, b), outcome(^, a, b)),
              "pow(", a, ", ", b, "): ", outcome(DM.det_pow, a, b), " vs Base ", outcome(^, a, b))
end
for (det_f, base_f) in ((DM.det_exp, exp), (DM.det_exp2, exp2), (DM.det_log, log), (DM.det_log2, log2), (DM.det_log10, log10)),
    x in (0.0, -0.0, 1.0, Inf, -Inf, NaN, -1.0, 1e-320, 800.0, -800.0)
    (d, b) = (outcome(det_f, x), outcome(base_f, x))
    @bp_check(isequal(d, b) || ((d isa Float64) && (b isa Float64) && (ulps(d, b) < 4)),
              det_f, "(", x, "): ", d, " vs Base ", b)
end
@bp_check(DM.det_pow(2, 10) === 1024) # Integer powers stay exact
@bp_check(DM.det_pow(0.5f0, 2.0f0) === 0.25f0)

# The real point: identical bits under different CPU targets, with and without FMA.
# The module has no dependencies, so it can be loaded standalone in a fresh process with any target.
if Sys.ARCH == :x86_64
    module_path = joinpath(@__DIR__, "..", "src", "DeterministicMath.jl")
    sweep_src = """
        include(raw"$module_path") # Also runs the self-check, which errors if the bits are wrong
        using .DeterministicMath
        h = UInt64(0)
        for i in 1:20_000
            t = i / 20_000
            for v in (det_exp(-700 + 1400t), det_exp2(-1070 + 2090t), det_log(t * 1e6), det_log2(t * 1e-300),
                      det_log10(t * 1e300), det_pow(t * 3.0, -20 + 40t), det_pow(-2.5 * t, Float64(i % 41 - 20)),
                      Float64(det_pow(Float32(t), 0.52f0)))
                global h = hash(reinterpret(UInt64, v), h)
            end
        end
        print(repr(h), " ", Core.Intrinsics.have_fma(Float64))
    """
    results = map(("generic", "sandybridge", "native")) do target
        cmd = `$(Base.julia_cmd()) --startup-file=no -C $target -e $sweep_src`
        target => read(pipeline(cmd; stdin=devnull), String)
    end
    hashes = Set(first(split(r)) for (_, r) in results)
    @bp_check(length(hashes) == 1, "DeterministicMath differs between CPU targets: ", results)
end

end)()
