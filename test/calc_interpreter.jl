# Tests for `CalcInterpreter`, the calculator language.

const CI = MJ.CalcInterpreter

(() -> begin

ENV_GRID = CI.CalcEnvironment(inputs = (:area, :scale), arrays = (:grid_size, ),
                              constants = (k = 3, ))

"Compiles and evaluates the given expression, checking that its source survives a round-trip"
function calc(expr, env = CI.CalcEnvironment(); rng = PRNG(0x1234), inputs...)
    program = CI.calc_compile(expr, env)
    @bp_check(CI.calc_compile(Meta.parse(string(program)), env) == program,
              "Round-trip failed: ", program, " vs ", Meta.parse(string(program)))
    return CI.calc_evaluate(program, rng; inputs...)
end
"Checks the value and its exact type (allowing for float error)"
function check_calc(expr, env::CI.CalcEnvironment, expected; kw...)
    actual = calc(expr, env; kw...)
    @bp_check((actual === expected) || (actual isa Float64 && expected isa Float64 && isapprox(actual, expected)),
              "`", expr, "` gave ", actual, " (", typeof(actual), "), expected ", expected, " (", typeof(expected), ")")
end
check_calc(expr, expected; kw...) = check_calc(expr, CI.CalcEnvironment(), expected; kw...)
"Checks that compiling or evaluating fails with a `CalcError` (but not its message, which may change)"
function check_calc_fails(expr, env = CI.CalcEnvironment(); inputs...)
    try
        CI.calc_evaluate(CI.calc_compile(expr, env), PRNG(0x1234); inputs...)
    catch e
        @bp_check(e isa CI.CalcError, "Expected a CalcError from `", expr, "`, got: ", sprint(showerror, e))
        return nothing
    end
    error("Should have failed: `", expr, "`")
end


# Literals and arithmetic, following Julia's types.
check_calc(:(1 + 2 * 3), 7)
check_calc(:(7 / 2), 3.5)
check_calc(:(7 ÷ 2), 3)
check_calc(:(sin(90)), 1.0) # Degrees
check_calc(:(cos(0)), 1.0)
check_calc(:(sinT(0.25) + cosT(0.5)), 0.0) # Full turns
# Bitwise operators, on Ints or Bools.
check_calc(:(6 & 3), 2)
check_calc(:(6 | 3), 7)
check_calc(:(1 | 2 | 4), 7)
check_calc(:(xor(6, 3)), 5)
check_calc(:(6 ⊻ 3), 5)
check_calc(:(~0), -1)
check_calc(:(1 << 4), 16)
check_calc(:(-16 >> 2), -4)
check_calc(:(-16 >>> 60), 15)
check_calc(:((5 & 4) == 4), true) # `&` binds tighter than comparisons
check_calc(:(true & false), false)
check_calc(:(true | false), true)
check_calc(:(~true), false)
check_calc_fails(:(1.5 & 1))
check_calc(:(7 % 3), 1)
check_calc(:(mod(-1, 3)), 2)
check_calc(:(2^10), 1024)
check_calc(:(2.0^-1), 0.5)
check_calc(:(-(1 + 2)), -3)
check_calc(:(2.5f0 * 2), 5.0)
check_calc(:(true + true), 2)
check_calc(5, 5)
check_calc(false, false)

# Constants, including user ones which replace built-ins.
check_calc(:(2pi), 2 * Float64(pi))
check_calc(:(π / 2), Float64(pi) / 2)
check_calc(:(k * 2), ENV_GRID, 6)
check_calc(:(pi), CI.CalcEnvironment(constants = Dict(:pi => 3)), 3)
check_calc(:(-Inf < 0), true)

# Built-in functions.
check_calc(:(clamp(15, 0, 10)), 10)
check_calc(:(sqrt(16)), 4.0)
check_calc(:(√16), 4.0)
check_calc(:(floor(2.7)), 2.0)
check_calc(:(ceil(2.1)), 3.0)
check_calc(:(floor(Int, 2.7)), 2)
check_calc(:(ceil(Int, 2.1)), 3)
check_calc(:(round(Int, 2.5)), 2)
check_calc(:(trunc(Int, -2.7)), -2)
check_calc(:(Int(4.0)), 4)
check_calc(:(float(4)), 4.0)
check_calc(:(min(5, 2, 8)), 2)
check_calc(:(max(5, 2.5)), 5.0)
check_calc(:(lerp(10, 20, 0.25)), 12.5)
check_calc(:(inv_lerp(10, 20, 12.5)), 0.25)
check_calc(:(exp(0)), 1.0)
check_calc(:(exp2(10)), 1024.0)
check_calc(:(log10(1000)), 3.0)
check_calc(:(log2(1 / 8)), -3.0)
check_calc(:(2.0 ^ 0.5), sqrt(2.0))
check_calc(:((-2.0) ^ 3), -8.0)
check_calc(:(2 ^ 10), 1024) # Int powers stay Int
check_calc(:(log(2, 8)), 3.0)
# Other trigonometry and friends aren't available (they'd need deterministic implementations).
for f in (:atan, :sinh, :hypot, :cbrt, :deg2rad)
    check_calc_fails(Expr(:call, f, 1))
end
check_calc_fails(:(log(-1)))
check_calc_fails(:((-8.0) ^ 0.5))
check_calc(:(iseven(4)), true)

# User functions, including one replacing a built-in.
env_funcs = CI.CalcEnvironment(functions = (double = CI.CalcFunction(x -> 2x, 1),
                                            clamp = CI.CalcFunction((x, a, b) -> -1, 3),
                                            anything = (args...) -> length(args)))
check_calc(:(double(4) + 1), env_funcs, 9)
check_calc(:(clamp(15, 0, 10)), env_funcs, -1)
check_calc(:(anything(1, 2, 3)), env_funcs, 3)
check_calc_fails(:(double(1, 2)), env_funcs)
check_calc_fails(:(clamp(1, 2)))
check_calc_fails(:(anything(1, 2, 3)))

# Comparisons and logic.
check_calc(:(1 < 2 <= 2 < 3), true)
check_calc(:(1 < 3 < 2), false)
check_calc(:(1 ≠ 2), true)
check_calc(:(!(1 == 1) || (2 ≥ 2 && 3 ≤ 4)), true)
check_calc(:(5 > 3 ? 10 : 20), 10)
check_calc(:(if 1 > 2
                 10
             elseif 2 > 1
                 20
             else
                 30
             end), 20)

# Inputs.
check_calc(:(area / 2 + scale), ENV_GRID, 51.5; area = 100, scale = 1.5)
check_calc(:(area / 2), ENV_GRID, 50.0; area = 100, grid_size = "ignored since it's unused")
check_calc_fails(:(area + scale), ENV_GRID; area = 100)
check_calc_fails(:(area), ENV_GRID; area = "big")

# Arrays and their generated functions, from tuples or B+ Vecs.
for grid_size in ((4, 9, 2, 9), Vec(4, 9, 2, 9))
    check_calc(:(grid_size(2)), ENV_GRID, 9; grid_size = grid_size)
    check_calc(:(grid_size_min()), ENV_GRID, 2; grid_size = grid_size)
    check_calc(:(grid_size_max()), ENV_GRID, 9; grid_size = grid_size)
    check_calc(:(grid_size_min_axis()), ENV_GRID, 3; grid_size = grid_size)
    check_calc(:(grid_size_max_axis()), ENV_GRID, 2; grid_size = grid_size)
    check_calc(:(grid_size_sum()), ENV_GRID, 24; grid_size = grid_size)
    check_calc(:(grid_size_prod()), ENV_GRID, 648; grid_size = grid_size)
    check_calc(:(grid_size_average()), ENV_GRID, 6.0; grid_size = grid_size)
    check_calc(:(grid_size_median()), ENV_GRID, 6.5; grid_size = grid_size)
    check_calc(:(grid_size_len()), ENV_GRID, 4; grid_size = grid_size)
    check_calc(:(grid_size(grid_size_len() - 1)), ENV_GRID, 2; grid_size = grid_size)
end
check_calc_fails(:(grid_size(5)), ENV_GRID; grid_size = (1, 2))
check_calc_fails(:(grid_size(1.0)), ENV_GRID; grid_size = (1, 2))
check_calc_fails(:(grid_size_min()), ENV_GRID; grid_size = ())
check_calc_fails(:(grid_size_min(1)), ENV_GRID)
check_calc_fails(:(grid_size), ENV_GRID)
check_calc_fails(:(grid_size(1)), ENV_GRID; grid_size = 5)

# Random intervals.
rng = PRNG(0xabcdef)
seen_ints = Set(calc(:(3:6); rng = rng) for _ in 1:200)
@bp_check(seen_ints == Set(3:6), seen_ints)
seen_reversed = Set(calc(:(6:3); rng = rng) for _ in 1:200)
@bp_check(seen_reversed == Set(3:6), seen_reversed)
floats = [ calc(:(2.0:0.5); rng = rng) for _ in 1:200 ]
@bp_check(all(f -> f isa Float64 && 0.5 <= f <= 2.0, floats), floats)
@bp_check(maximum(floats) - minimum(floats) > 1.0, "Floats should spread across the interval")
bools = Set(calc(:(true:false); rng = rng) for _ in 1:100)
@bp_check(bools == Set([ false, true ]), bools)
@bp_check(all(_ -> calc(:(true:true); rng = rng) === true, 1:20))
check_calc(:(5:5), 5)
check_calc(:(2 * (1:1) + 1), 3)
check_calc(:(area:area), ENV_GRID, 7; area = 7)
# Same seed, same results.
@bp_check([ calc(:(1:1000); rng = PRNG(0x55)) for _ in 1:1 ] == [ calc(:(1:1000); rng = PRNG(0x55)) ])
# Mixing Int and Float64 makes a Float64 interval.
mixed = [ calc(:(1:2.5); rng = rng) for _ in 1:200 ]
@bp_check(all(f -> f isa Float64 && 1.0 <= f <= 2.5, mixed), mixed)
@bp_check(length(unique(mixed)) > 100, "Mixed intervals should be continuous")
@bp_check(calc(:(1:area), ENV_GRID; rng = rng, area = 2.5) isa Float64)
@bp_check(calc(:(area:3), ENV_GRID; rng = rng, area = 3.0) === 3.0)
@bp_check(calc(:(curve(1:2.0, 2)); rng = rng) isa Float64)
# Mixing Bool with numbers fails while compiling if the ends are constant, otherwise while evaluating.
check_calc_fails(:(false:1))
check_calc_fails(:(1.5:true))
check_calc_fails(:(1:area), ENV_GRID; area = true)
check_calc_fails(:(curve(area:1, 2)), ENV_GRID; area = false)
check_calc_fails(:(1:2:5))

# Curved random intervals lean toward the first end when the power is above 1, and the second when below.
mean_of(src; n = 4000, kw...) = let rng = PRNG(0x9876)
    sum(_ -> calc(src, ENV_GRID; rng = rng, kw...), 1:n) / n
end
@bp_check(isapprox(mean_of(:(curve(0.0:1.0, 2))), 1/3; atol=0.03), mean_of(:(curve(0.0:1.0, 2))))
@bp_check(isapprox(mean_of(:(curve(1.0:0.0, 2))), 2/3; atol=0.03), mean_of(:(curve(1.0:0.0, 2))))
@bp_check(isapprox(mean_of(:(curve(0.0:1.0, 0.5))), 2/3; atol=0.03), mean_of(:(curve(0.0:1.0, 0.5))))
@bp_check(isapprox(mean_of(:(curve(0.0:1.0, 1))), 1/2; atol=0.03), mean_of(:(curve(0.0:1.0, 1))))
# Ints use equal bins from the first end, so a power of 1 is uniform.
ints = [ calc(:(curve(1:10, 3)); rng = rng) for _ in 1:2000 ]
@bp_check(all(i -> i isa Int && 1 <= i <= 10, ints) && Set(ints) == Set(1:10), Set(ints))
@bp_check(sum(ints) / length(ints) < 3.5, sum(ints) / length(ints))
ints = [ calc(:(curve(10:1, 3)); rng = rng) for _ in 1:2000 ]
@bp_check(Set(ints) == Set(1:10) && sum(ints) / length(ints) > 7.5, sum(ints) / length(ints))
@bp_check(isapprox(mean_of(:(curve(1:10, 1))), 5.5; atol=0.2), mean_of(:(curve(1:10, 1))))
check_calc(:(curve(5:5, 2)), 5)
check_calc(:(curve(true:true, 0.1)), true)
# Bools are the two-bin case: `t^2 >= 0.5` happens with probability 1 - sqrt(0.5).
@bp_check(isapprox(mean_of(:(curve(false:true, 2))), 1 - sqrt(0.5); atol=0.03), mean_of(:(curve(false:true, 2))))
# The power can come from inputs.
@bp_check(isapprox(mean_of(:(curve(0.0:1.0, area)); area = 3), 1/4; atol=0.03))
# Errors.
check_calc_fails(:(curve(1, 2)))
check_calc_fails(:(curve(1:2)))
check_calc_fails(:(curve(1:false, 2)))
check_calc_fails(:(curve(1:2, 0)))
check_calc_fails(:(curve(1:2, true)))
check_calc_fails(:(curve(1:2, -area)), ENV_GRID; area = 1)
@bp_check(calc(:(curve(1:area, 2)), ENV_GRID; area = 1.5) isa Float64)

# Constant folding.
@bp_check(CI.calc_is_constant(CI.calc_compile(:(2pi + k * 2), ENV_GRID)))
@bp_check(CI.calc_constant_value(CI.calc_compile(:(k * 2), ENV_GRID)) === 6)
@bp_check(!CI.calc_is_constant(CI.calc_compile(:(1:2))))
@bp_check(!CI.calc_is_constant(CI.calc_compile(:(area * 2), ENV_GRID)))
@bp_check(CI.calc_constant_value(CI.calc_compile(:(1 > 0 ? 3 : area), ENV_GRID)) === 3)
# Errors in constant branches that are never taken don't matter...
check_calc(:(false ? div(1, 0) : 3), 3)
check_calc(:(area > 0 ? 3 : div(1, 0)), ENV_GRID, 3; area = 1)
# ...but they do when evaluated.
check_calc_fails(:(div(1, 0)))
check_calc_fails(:(sqrt(-1)))
check_calc_fails(:(area > 0 ? div(1, 0) : 3), ENV_GRID; area = 1)

# Runtime errors can be handled with a fallback value, and name the calculation.
handled_msgs = String[ ]
log_and_fallback(e) = (push!(handled_msgs, e.msg); -1)
@bp_check(CI.calc_evaluate(CI.calc_compile(:(div(area, 0)), ENV_GRID), PRNG(1), log_and_fallback; area = 3) === -1)
@bp_check(CI.calc_evaluate(CI.calc_compile(:(area + 1), ENV_GRID), PRNG(1), log_and_fallback) === -1)
@bp_check(length(handled_msgs) == 2, handled_msgs)
# Errors raised by user functions are handled too.
env_erroring = CI.CalcEnvironment(inputs = (:x, ), functions = (fussy = CI.CalcFunction(v -> error("too fussy"), 1), ))
@bp_check(CI.calc_evaluate(CI.calc_compile(:(fussy(x)), env_erroring), PRNG(1), e -> 2.5; x = 1) === 2.5)
check_calc_fails(:(fussy(x)), env_erroring; x = 1)
# The handler's own result must be a number.
try
    CI.calc_evaluate(CI.calc_compile(:(div(1, area)), ENV_GRID), PRNG(1), e -> "oops"; area = 0)
    error("Should have failed")
catch e
    @bp_check(e isa CI.CalcError, sprint(showerror, e))
end

# Type errors.
check_calc_fails(:(1 ? 2 : 3))
check_calc_fails(:(1 && true))
check_calc_fails(:(area ? 2 : 3), ENV_GRID; area = 1)
check_calc_fails(:(bad()),
                 CI.CalcEnvironment(functions = (bad = CI.CalcFunction(() -> "hi", 0), )))

# Unsupported syntax and unknown names.
check_calc_fails(:(x + 1))
check_calc_fails(:(foo(1)))
check_calc_fails(:(x = 1))
check_calc_fails(:(a.b))
check_calc_fails(:([1, 2]))
check_calc_fails(:((1, 2)))
check_calc_fails(:(max(1; init=2)))
check_calc_fails(:(max(xs...)))
check_calc_fails(:(begin 1; 2 end))
check_calc_fails(:(if true 1 end))
check_calc_fails("hello")
check_calc_fails(:(nothing))
check_calc_fails(:(floor(Float64, 2.5)))

# Name clashes in the environment.
for kw in (
    (constants = (a = 1, ), inputs = (:a, )),
    (arrays = (:min, ), ),
    (arrays = (:grid, ), functions = (grid_sum = sum, )),
    (arrays = (:grid, :grid_min), ),
    (functions = (curve = sum, ), ),
    (arrays = (:curve, ), ),
)
    try
        CI.CalcEnvironment(; kw...)
        error("Should have failed: ", kw)
    catch e
        @bp_check(e isa CI.CalcError, sprint(showerror, e))
    end
end

# Equality and printing go by the source.
@bp_check(CI.calc_compile(:(1 + area), ENV_GRID) == CI.calc_compile(:(1 + area), ENV_GRID))
@bp_check(CI.calc_compile(:(1 + area), ENV_GRID) != CI.calc_compile(:(area + 1), ENV_GRID))
@bp_check(string(CI.calc_compile(:(1 + area), ENV_GRID)) == "1 + area")


#################
# User functions

ENV_MACROS = CI.CalcEnvironment(inputs = (:area, :scale), arrays = (:grid_size, ),
                                constants = (k = 3, ), macros = true)
macros_of(defs...) = Dict{Symbol, CI.CalcMacro}(CI.calc_macro_parse.(defs)...)
"Compiles in `ENV_MACROS` and evaluates with the given user functions"
eval_macros(expr, macros; rng = PRNG(0x1234), inputs...) =
    CI.calc_evaluate(CI.calc_compile(expr, ENV_MACROS), rng, CI.calc_throw, macros; inputs...)
function check_macros_fail(expr, macros; inputs...)
    try
        eval_macros(expr, macros; inputs...)
    catch e
        @bp_check(e isa CI.CalcError, "Expected a CalcError from `", expr, "`, got: ", sprint(showerror, e))
        return nothing
    end
    error("Should have failed: `", expr, "`")
end

# Parsing a definition, from an expression or a string, and printing it back.
let (name, m) = CI.calc_macro_parse(:( half(a) = a / 2 ))
    @bp_check(name == :half)
    @bp_check(m.params == [ :a ])
    @bp_check(m.body == :(a / 2), m.body) # Line numbers are stripped
    @bp_check(CI.calc_macro_parse("half(a) = a / 2") == (name => m))
    @bp_check(CI.calc_macro_parse(CI.calc_macro_expr(name, m)) == (name => m))
end
@bp_check(CI.calc_macro_parse("function two() 2 end") == (:two => CI.CalcMacro(Symbol[ ], 2)))

# Definitions only check their shape: unknown names in the body are fine until it runs.
CI.calc_macro_parse(:( later() = no_such_thing + also_missing(1) ))
for def in (:( 5 ), :( f ), :( f = 3 ), :( a.b(x) = x ), :( f(1) = 2 ), :( f(a, a) = a ),
            :( clamp(a) = a ), :( curve(a) = a ))
    try
        CI.calc_macro_parse(def)
        error("Should have failed: ", def)
    catch e
        @bp_check(e isa CI.CalcError, "Expected a CalcError from `", def, "`, got: ", sprint(showerror, e))
    end
end

# Calls are resolved against the table given when evaluating, and their bodies use the caller's environment.
MACROS = macros_of(:( half(a) = a / 2 ), :( cells() = grid_size_prod() ), :( kk() = k * k ),
                   :( pick(c, a, b) = c ? a : b ), :( both(a, b) = a + b ))
@bp_check(eval_macros(:( half(10) ), MACROS) === 5.0)
@bp_check(eval_macros(:( half(area) + 1 ), MACROS; area = 4) === 3.0)
@bp_check(eval_macros(:( cells() ), MACROS; grid_size = (3, 4)) === 12) # Input only the body uses
@bp_check(eval_macros(:( kk() ), MACROS) === 9) # Environment constant
@bp_check(eval_macros(:( pick(area > 2, 10, 20) ), MACROS; area = 3) === 10)
@bp_check(eval_macros(:( half(half(8)) ), MACROS) === 2.0) # Nesting
@bp_check(eval_macros(:( both(half(2), kk()) ), MACROS) === 10.0)

# Each argument is evaluated once, by the caller, so a random argument is the same everywhere it's used.
let m = macros_of(:( twice(a) = a + a ), :( twice_minus(a) = a - a ))
    results = Set(eval_macros(:( twice(1:6) ), m; rng = PRNG(i)) for i in 1:60)
    @bp_check(results == Set(2:2:12), results)
    @bp_check(all(i -> eval_macros(:( twice_minus(1:1000) ), m; rng = PRNG(i)) == 0, 1:20))
    # Arguments are evaluated even if the body never uses them.
    @bp_check(eval_macros(:( pick(true, 1, 1:6) ), macros_of(:( pick(c, a, b) = c ? a : b ))) === 1)
end
# Parameters hide the environment's constants and inputs with the same name.
let m = macros_of(:( shadow_k(k) = k + 1 ), :( shadow_area(area) = area * 2 ))
    @bp_check(eval_macros(:( shadow_k(10) ), m) === 11)
    @bp_check(eval_macros(:( shadow_area(1) + area ), m; area = 4) === 6)
end
# A function can't see the parameters of the function calling it.
check_macros_fail(:( outer(1) ), macros_of(:( outer(p) = inner() ), :( inner() = p )))
# Parameters which share a name with something else work too.
let m = macros_of(:( f(scale, max) = max(scale, 1) ), :( g(a) = 0 < a <= 2 ), :( h(a) = f(a, 0) ))
    @bp_check(eval_macros(:( f(5, 0) ), m) === 5) # Call heads aren't parameters
    @bp_check(eval_macros(:( g(2) ), m) === true)
    @bp_check(eval_macros(:( h(scale) ), m; scale = 7) === 7) # `scale` passes through `f`'s parameter of the same name
end

# Finding out which inputs a calculation uses, even in branches that are never taken or in user functions.
let m = macros_of(:( uses_area() = area + 1 ), :( shadows_area(area) = area ), :( indirect() = uses_area() ),
                  :( loops(a) = loops(a) ))
    mentions_area(expr) = CI.calc_mentions(CI.calc_compile(expr, ENV_MACROS), (:area, ), m)
    @bp_check(mentions_area(:( area * 2 )))
    @bp_check(!mentions_area(:( scale * 2 )))
    @bp_check(mentions_area(:( false ? area : 1 ))) # Folded away while compiling, but still mentioned
    @bp_check(mentions_area(:( 1 < area <= 3 )))
    @bp_check(mentions_area(:( indirect() )))
    @bp_check(!mentions_area(:( shadows_area(3) )))  # The parameter hides the input
    @bp_check(mentions_area(:( shadows_area(area) ))) # ...but the argument is the caller's
    @bp_check(!mentions_area(:( loops(1) )))
    @bp_check(!mentions_area(:( not_defined() )))
end

# The table may change between evaluations without recompiling the caller.
let program = CI.calc_compile(:( 1 + bonus() ), ENV_MACROS),
    table = macros_of("bonus() = 1")
    @bp_check(!CI.calc_is_constant(program))
    @bp_check(CI.calc_evaluate(program, PRNG(1), CI.calc_throw, table) === 2)
    push!(table, CI.calc_macro_parse("bonus() = 10"))
    @bp_check(CI.calc_evaluate(program, PRNG(1), CI.calc_throw, table) === 11)
    delete!(table, :bonus)
    @bp_check(CI.calc_evaluate(program, PRNG(1), e -> -1, table) === -1)
end

# The caller's arguments are checked when compiling; everything else waits until the call runs.
check_calc_fails(:( later(missing_name) ), ENV_MACROS)
CI.calc_compile(:( later(1) + also_later() ), ENV_MACROS)
check_macros_fail(:( nope(1) ), MACROS)
check_macros_fail(:( half(1, 2) ), MACROS)
check_macros_fail(:( broken() ), macros_of("broken() = missing_name"))
check_macros_fail(:( wants_len() ), macros_of("wants_len() = len")) # Not in this environment
check_macros_fail(:( bad_axis() ), macros_of("bad_axis() = grid_size(5)");
                  grid_size = (1, 2))
check_macros_fail(:( forever(1) ), macros_of("forever(a) = forever(a + 1)"))
# Without `macros = true`, unknown functions stay a compile error.
check_calc_fails(:( half(1) ), ENV_GRID)

println("calc_interpreter tests passed")

end)()
