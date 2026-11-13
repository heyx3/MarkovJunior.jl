"
A small calculator language: a subset of Julia expressions,
  checked and compiled once, then evaluated any number of times without invoking Julia's compiler.
It's a plain number-cruncher; everything specific to its user is injected through a `CalcEnvironment`.

(maintained mostly by Claude because I didn't want to write a whole calculator parser+evaluator)

Values are `Bool`, `Int`, or `Float64` (see `CalcValue`), with Julia's own arithmetic rules
  (e.g. `7/2 == 3.5` while `7÷2 == 3`).

Supported syntax:
* Number and `Bool` literals
* Named constants, either built-in (like `pi`) or injected by the user
* Named inputs, whose values are given each time the calculation is evaluated
* Function calls, either built-in (like `clamp` and `sqrt`) or injected by the user.
  `floor`, `ceil`, `round`, and `trunc` also accept `Int` as their first argument, like in Julia.
* Named 1D arrays given each time the calculation is evaluated.
  An array called `a` offers the functions `a(i)`, `a_min()`, `a_max()`, `a_min_axis()`, `a_max_axis()`,
    `a_sum()`, `a_prod()`, `a_average()`, `a_median()`, and `a_len()`.
* Comparisons, including chains like `0 < x <= 1`
* `!`, `&&`, `||`, `a ? b : c`, and `if`/`elseif`/`else` (which must have an `else`)
* Random intervals `a:b`, which pick a random value from `a` to `b` (inclusive) every time they're evaluated.
  Two `Int`s pick an `Int`, while any other pair of numbers picks a `Float64` (so `0.5:2` is a `Float64` interval).
  `false:true` flips a coin, and `Bool`s can't be mixed with numbers. The ends may come in either order.
* Curved random intervals `curve(a:b, p)`, which lean toward `a` when `p > 1` and toward `b` when `p < 1`.
  A uniform `t` from 0 to 1 is raised to the power `p`, then mapped from `a` to `b`
    (for `Int` and `Bool` intervals, into equal-sized bins starting from `a`).
  The first argument must be written as an interval; `curve` is reserved for this and isn't a normal function.
* User functions (`CalcMacro`), if the environment allows them: calls to unknown functions are looked up
    by name each time they're evaluated, in a table given to `calc_evaluate()`.

Compile with `calc_compile(expr, env)`, then evaluate with `calc_evaluate(program, rng[, on_error[, macros]]; inputs...)`.

Printing a `CalcProgram` gives source which behaves the same, but isn't necessarily the original text.
"
module CalcInterpreter

using Random, Statistics, PrecompileTools
# Results must be identical on every machine, so math functions come from here rather than Base.
using ..DeterministicMath

export CalcValue, CalcError,
       CalcFunction, CalcEnvironment, CalcProgram,
       CalcMacro, CalcMacros, calc_macro_parse, calc_macro_expr,
       calc_compile, calc_evaluate, calc_throw,
       calc_is_constant, calc_constant_value, calc_mentions


#########################
#  Values and errors

"The types of value that a calculation works with"
const CalcValue = Union{Bool, Int, Float64}

"An error in a calculation, either while compiling it or while evaluating it"
struct CalcError <: Exception
    msg::String
end
Base.showerror(io::IO, e::CalcError) = print(io, "CalcError: ", e.msg)

calc_error(msg...) = throw(CalcError(string(msg...)))

"Converts a Julia number into a `CalcValue`, or throws a `CalcError`"
to_calc_value(x::Bool, context = nothing)::CalcValue = x
to_calc_value(x::Integer, context = nothing)::CalcValue = if typemin(Int) <= x <= typemax(Int)
    convert(Int, x)
else
    calc_error("Integer `", x, "` is too large", context_suffix(context))
end
to_calc_value(x::Union{AbstractFloat, AbstractIrrational}, context = nothing)::CalcValue = convert(Float64, x)
to_calc_value(x, context = nothing)::CalcValue = calc_error(
    "Calculations only work with Bool, Int, and Float64 values; got a ", typeof(x), " (`", x, "`)",
    context_suffix(context)
)

context_suffix(::Nothing) = ""
context_suffix(context) = string(", from ", context)

"Converts a Julia array or tuple into a vector of `CalcValue`, or throws a `CalcError`"
to_calc_array(x::Union{AbstractVector, Tuple}, name::Symbol)::Vector{CalcValue} =
    CalcValue[ to_calc_value(e, "array `$name`") for e in x ]
to_calc_array(x, name::Symbol) = calc_error("Input `", name, "` should be a 1D array or tuple; got a ", typeof(x))


###############
#  Functions

const VARIADIC = typemax(Int)

"
A Julia function made available to calculations.
`n_args` is the allowed range for the number of arguments, checked while compiling.

The function should be pure (always giving the same output for the same inputs),
  because calls with constant arguments are evaluated once while compiling.
It should return a number (see `CalcValue`), and is allowed to throw common math errors
  like `DomainError` or `DivideError` (or call `error()`), which become a `CalcError`.
"
struct CalcFunction
    f::Any
    n_args::UnitRange{Int}
end
CalcFunction(f, n_args::Integer) = CalcFunction(f, n_args:n_args)
CalcFunction(f, min_args::Integer, max_args::Integer) = CalcFunction(f, min_args:max_args)
"Without an argument count, the function accepts any number of arguments (checked only when it's called)"
CalcFunction(f) = CalcFunction(f, 0:VARIADIC)

format_arg_count(r::UnitRange{Int}) = if first(r) == last(r)
    string(first(r))
elseif last(r) == VARIADIC
    string(first(r), " or more")
else
    string(first(r), " to ", last(r))
end

"Calls a function from inside a calculation, turning its common errors into `CalcError`"
function call_calc_function(f::CalcFunction, name::Symbol, args::CalcValue...)::CalcValue
    result = try
        f.f(args...)
    catch e
        if e isa Union{DomainError, DivideError, InexactError, OverflowError,
                       ArgumentError, BoundsError, MethodError, ErrorException}
            calc_error("`", name, "(", join(args, ", "), ")` failed: ", sprint(showerror, e))
        else
            rethrow()
        end
    end
    return to_calc_value(result, "function `$name`")
end

lerp(a, b, t) = a + ((b - a) * t)
inv_lerp(a, b, x) = (x - a) / (b - a)
"`log(x)`, or `log(b, x)` for a logarithm with base `b`, like Base"
det_log_any_base(x) = det_log(x)
det_log_any_base(b, x) = det_log(x) / det_log(b)

const BUILTIN_FUNCTIONS = Dict{Symbol, CalcFunction}(
    # Arithmetic:
    :+ => CalcFunction(+, 1, VARIADIC),
    :- => CalcFunction(-, 1, 2),
    :* => CalcFunction(*, 2, VARIADIC),
    :/ => CalcFunction(/, 2),
    :^ => CalcFunction(det_pow, 2),
    :÷ => CalcFunction(div, 2),
    :% => CalcFunction(rem, 2),
    :div => CalcFunction(div, 2),
    :rem => CalcFunction(rem, 2),
    :mod => CalcFunction(mod, 2),
    :fld => CalcFunction(fld, 2),
    :cld => CalcFunction(cld, 2),
    :abs => CalcFunction(abs, 1),
    :sign => CalcFunction(sign, 1),
    :min => CalcFunction(min, 1, VARIADIC),
    :max => CalcFunction(max, 1, VARIADIC),
    :clamp => CalcFunction(clamp, 3),
    :lerp => CalcFunction(lerp, 3),
    :inv_lerp => CalcFunction(inv_lerp, 3),

    # Rounding and conversion:
    :floor => CalcFunction(floor, 1),
    :ceil => CalcFunction(ceil, 1),
    :round => CalcFunction(round, 1),
    :trunc => CalcFunction(trunc, 1),
    :Int => CalcFunction(Int, 1),
    :Float64 => CalcFunction(Float64, 1),
    :float => CalcFunction(float, 1),
    :iseven => CalcFunction(iseven, 1),
    :isodd => CalcFunction(isodd, 1),

    # Powers, logarithms, and trigonometry.
    # These come from `DeterministicMath`, as Base's versions can differ between machines in the last bit.
    # `sqrt` is exactly rounded by IEEE 754, so Base's is already deterministic.
    :sqrt => CalcFunction(sqrt, 1),
    :√ => CalcFunction(sqrt, 1),
    :exp => CalcFunction(det_exp, 1),
    :exp2 => CalcFunction(det_exp2, 1),
    :log => CalcFunction(det_log_any_base, 1, 2),
    :log2 => CalcFunction(det_log2, 1),
    :log10 => CalcFunction(det_log10, 1),
    # Angles are in degrees, or in full turns for the `T` versions (which have a period of 1, like a `lerp` `t`).
    :sin => CalcFunction(det_sind, 1),
    :cos => CalcFunction(det_cosd, 1),
    :tan => CalcFunction(det_tand, 1),
    :sinT => CalcFunction(det_sin_turns, 1),
    :cosT => CalcFunction(det_cos_turns, 1),

    # Bitwise, on Ints (or on Bools, as logic which doesn't short-circuit):
    :& => CalcFunction(&, 2),
    :| => CalcFunction(|, 2),
    :xor => CalcFunction(xor, 2),
    :⊻ => CalcFunction(xor, 2),
    :~ => CalcFunction(~, 1),
    :<< => CalcFunction(<<, 2),
    :>> => CalcFunction(>>, 2),
    :>>> => CalcFunction(>>>, 2),

    # Logic and comparison:
    :! => CalcFunction(!, 1),
    :(==) => CalcFunction(==, 2),
    :!= => CalcFunction(!=, 2),
    :≠ => CalcFunction(!=, 2),
    :< => CalcFunction(<, 2),
    :<= => CalcFunction(<=, 2),
    :≤ => CalcFunction(<=, 2),
    :> => CalcFunction(>, 2),
    :>= => CalcFunction(>=, 2),
    :≥ => CalcFunction(>=, 2)
)

"Versions of the rounding functions which take `Int` as their first argument, like `floor(Int, x)`"
const BUILTIN_ROUNDING_TO_INT = Dict{Symbol, CalcFunction}(
    :floor => CalcFunction(x -> floor(Int, x), 1),
    :ceil => CalcFunction(x -> ceil(Int, x), 1),
    :round => CalcFunction(x -> round(Int, x), 1),
    :trunc => CalcFunction(x -> trunc(Int, x), 1)
)

"Function-call syntax with special meaning, which users can't define"
const RESERVED_FUNCTION_NAMES = Set{Symbol}([ :curve ])

const BUILTIN_CONSTANTS = Dict{Symbol, CalcValue}(
    :pi => Float64(pi),
    :π => Float64(pi),
    :ℯ => Float64(ℯ),
    :Inf => Inf
)

"
The functions offered for each named array, by their suffix.
The array `a` also offers `a(i)` to get an element.
"
const ARRAY_QUERIES = Dict{String, Any}(
    "_min" => minimum,
    "_max" => maximum,
    "_min_axis" => argmin,
    "_max_axis" => argmax,
    "_sum" => sequential_sum,
    "_prod" => sequential_prod,
    "_average" => sequential_mean,
    "_median" => median,
    "_len" => length
)


############
#  Macros

"
A user-defined function, written in the calculator language like `half(a) = a / 2`
  (constants are just functions without parameters, like `k() = 3`).

Its arguments are evaluated once, by the caller, like a normal function's
  (so with `twice(a) = a + a`, `twice(1:6)` is always even).
But its body is like a macro's:
* Defining one needs no `CalcEnvironment`, and its body isn't checked when it's defined.
  It's compiled when a call to it is first evaluated, in the environment of the calculation making that call
    (so the body may use that environment's inputs, arrays, and functions, as well as its own parameters).
* Calls look it up by name every time they're evaluated, in the table passed to `calc_evaluate()`,
    so macros may be added or replaced between evaluations without recompiling their callers.

Create one from its definition with `calc_macro_parse()`.
"
struct CalcMacro
    params::Vector{Symbol}
    body::Any # Normalized like a calculation's expression
end
Base.:(==)(a::CalcMacro, b::CalcMacro) = (a.params == b.params) && (a.body == b.body)
Base.hash(m::CalcMacro, h::UInt) = hash(m.body, hash(m.params, hash(CalcMacro, h)))

"A table of user functions, by name"
const CalcMacros = AbstractDict{Symbol, CalcMacro}
const NO_CALC_MACROS = Dict{Symbol, CalcMacro}()

"Calls to user functions can't nest deeper than this, so infinite recursion becomes a `CalcError`"
const MAX_MACRO_DEPTH = 64

"
Parses a user function's definition, like `:( half(a) = a / 2 )` or the string `\"half(a) = a / 2\"`,
  returning `name => CalcMacro`.
Only the definition's shape is checked here; the body is checked when a call to it runs.
Throws a `CalcError` if the definition is malformed.
"
calc_macro_parse(src::AbstractString)::Pair{Symbol, CalcMacro} = calc_macro_parse(Meta.parse(src))
function calc_macro_parse(def)::Pair{Symbol, CalcMacro}
    if !Base.isexpr(def, (:(=), :function), 2) || !Base.isexpr(def.args[1], :call)
        calc_error("A user function is defined like `name(a, b) = a + b`; got `", def, "`")
    end
    (name, params...) = def.args[1].args
    if !isa(name, Symbol)
        calc_error("A user function's name must be a plain name; got `", name, "` in `", def, "`")
    elseif haskey(BUILTIN_FUNCTIONS, name) || (name in RESERVED_FUNCTION_NAMES) || (name == :(:))
        calc_error("`", name, "` is built into the calculator language, so it can't be a user function")
    end
    for p in params
        if !isa(p, Symbol)
            calc_error("User function `", name, "` has a parameter which isn't a plain name: `", p, "`")
        elseif count(==(p), params) > 1
            calc_error("User function `", name, "` has the parameter `", p, "` more than once")
        end
    end
    return name => CalcMacro(collect(Symbol, params), normalize_expr(def.args[2]))
end

"Gets the source of a user function's definition, which `calc_macro_parse()` turns back into the same macro"
calc_macro_expr(name::Symbol, m::CalcMacro)::Expr = Expr(:(=), Expr(:call, name, m.params...), m.body)


#################
#  Environment

"
Everything a calculation can reference, on top of the language's built-ins:

* `constants`: named values which are known while compiling (e.g. `(k = 3, )`).
They may replace built-in constants like `pi`.
* `functions`: named `CalcFunction`s (or plain Julia functions with any number of arguments).
They may replace built-in functions like `clamp`.
* `inputs`: names of the values which are given each time the calculation is evaluated.
* `arrays`: names of the 1D arrays which are given each time the calculation is evaluated.
Each array also offers some generated functions (see `ARRAY_QUERIES`), whose names must not clash with anything.

`constants` and `functions` can be given as a `NamedTuple`, a `Dict`, or a collection of `Pair`s.

If `macros` is true, calling an unknown function isn't an error while compiling:
  it's assumed to be a user function (see `CalcMacro`), looked up when the call is evaluated.
Otherwise it's a compile error, as the environment's user doesn't offer any user functions.
"
struct CalcEnvironment
    constants::Dict{Symbol, CalcValue}
    functions::Dict{Symbol, CalcFunction}
    inputs::Set{Symbol}
    arrays::Set{Symbol}
    # Each generated array function, mapped to its array and query (`nothing` for element access).
    array_functions::Dict{Symbol, Tuple{Symbol, Any}}
    macros::Bool
end
"
The environment a user function's body is compiled in: the caller's environment,
  plus the function's parameters as inputs (which hide any constants or inputs with the same names)
"
calc_macro_body_env(env::CalcEnvironment, params::Vector{Symbol}) = CalcEnvironment(
    filter(kvp -> !in(kvp[1], params), env.constants),
    env.functions,
    union(env.inputs, params),
    env.arrays,
    env.array_functions,
    env.macros
)
function CalcEnvironment(; constants = (), functions = (), inputs = (), arrays = (), macros::Bool = false)
    all_constants = copy(BUILTIN_CONSTANTS)
    for (name, value) in named_pairs(constants)
        all_constants[Symbol(name)] = to_calc_value(value, "constant `$name`")
    end

    all_functions = copy(BUILTIN_FUNCTIONS)
    for (name, f) in named_pairs(functions)
        if Symbol(name) in RESERVED_FUNCTION_NAMES
            calc_error("`", name, "` is reserved by the calculator language and can't be a user function")
        end
        all_functions[Symbol(name)] = (f isa CalcFunction) ? f : CalcFunction(f)
    end

    input_set = Set{Symbol}(Symbol.(inputs))
    for name in input_set
        if haskey(all_constants, name)
            calc_error("`", name, "` can't be both a constant and an input")
        end
    end

    array_set = Set{Symbol}(Symbol.(arrays))
    array_functions = Dict{Symbol, Tuple{Symbol, Any}}()
    for array_name in array_set
        for (fn_name, query) in Iterators.flatten((
            ((array_name, nothing), ),
            ((Symbol(array_name, suffix), q) for (suffix, q) in ARRAY_QUERIES)
        ))
            if haskey(all_functions, fn_name) || (fn_name in RESERVED_FUNCTION_NAMES)
                calc_error("Array `", array_name, "` offers a function `", fn_name, "`, ",
                             "which clashes with an existing function")
            elseif haskey(array_functions, fn_name)
                calc_error("Array `", array_name, "` offers a function `", fn_name, "`, ",
                             "which clashes with one from array `", array_functions[fn_name][1], "`")
            end
            array_functions[fn_name] = (array_name, query)
        end
    end

    return CalcEnvironment(all_constants, all_functions, input_set, array_set, array_functions, macros)
end

"Iterates over `name => value` pairs, from a `NamedTuple`, a `Dict`, or a collection of `Pair`s"
named_pairs(x::Union{NamedTuple, AbstractDict}) = pairs(x)
named_pairs(x) = x


###########
#  Nodes

abstract type CalcNode end

struct CalcNodeConst <: CalcNode
    value::CalcValue
end
struct CalcNodeInput <: CalcNode
    slot::Int
end
struct CalcNodeArrayElement <: CalcNode
    name::Symbol
    slot::Int
    axis::CalcNode
end
struct CalcNodeArrayQuery <: CalcNode
    name::Symbol # The generated function's name
    slot::Int
    query::Any
end
struct CalcNodeCall <: CalcNode
    name::Symbol
    f::CalcFunction
    args::Vector{CalcNode}
end
"A chain of comparisons like `a < b <= c`, which stops at the first `false`"
struct CalcNodeComparison <: CalcNode
    args::Vector{CalcNode}
    op_names::Vector{Symbol}
    ops::Vector{CalcFunction}
end
struct CalcNodeAnd <: CalcNode
    a::CalcNode
    b::CalcNode
end
struct CalcNodeOr <: CalcNode
    a::CalcNode
    b::CalcNode
end
struct CalcNodeIf <: CalcNode
    condition::CalcNode
    if_true::CalcNode
    if_false::CalcNode
end
struct CalcNodeInterval <: CalcNode
    a::CalcNode
    b::CalcNode
end
"A random interval whose uniform sample is raised to some power, leaning toward one end"
struct CalcNodeCurvedInterval <: CalcNode
    a::CalcNode
    b::CalcNode
    power::CalcNode
end
"
A call to a user function (see `CalcMacro`), which is looked up by name every time it's evaluated.
Its body is compiled in the caller's environment plus its parameters, and cached until the macro is replaced.
"
struct CalcNodeMacroCall <: CalcNode
    name::Symbol
    args::Vector{CalcNode} # Evaluated by the caller, then given to the body as its parameters
    arg_exprs::Vector{Any} # For error messages
    # The environment of the outermost calculation, without any user functions' parameters,
    #   so a user function never sees the parameters of the one calling it.
    env::CalcEnvironment
    # The last compiled body, as `(macro, program)`.
    # Kept as one immutable tuple so concurrent evaluations never see a mismatched pair.
    cache::Base.RefValue{Any}
end


#############
#  Program

"A compiled calculation, ready to evaluate with `calc_evaluate()`"
struct CalcProgram
    # The original expression, minus line numbers.
    expr::Any
    root::CalcNode
    # The inputs and arrays which the calculation actually uses, in the order of their slots.
    input_names::Vector{Symbol}
    array_names::Vector{Symbol}
end

# Programs are compared and printed by their source;
#    note that the same source compiled in different environments compares equal.
Base.:(==)(a::CalcProgram, b::CalcProgram) = (a.expr == b.expr)
Base.hash(p::CalcProgram, h::UInt) = hash(p.expr, hash(CalcProgram, h))
Base.print(io::IO, p::CalcProgram) = print(io, p.expr)
Base.show(io::IO, p::CalcProgram) = print(io, "CalcProgram(`", p.expr, "`)")

"Whether the calculation always gives the same value (which is then available from `calc_constant_value()`)"
calc_is_constant(p::CalcProgram)::Bool = p.root isa CalcNodeConst
calc_constant_value(p::CalcProgram)::CalcValue = if calc_is_constant(p)
    p.root.value
else
    error("Calculation `", p.expr, "` isn't constant")
end

"
Finds whether the calculation mentions any of the given input (or array) names,
  either directly or in the body of any nested user function calls.

Works by examining the source, so it intentionally includes dead code
  which is normally optimized away.
"
function calc_mentions(p::CalcProgram, names, macros::CalcMacros = NO_CALC_MACROS)::Bool
    checked_macros = Set{Symbol}()
    mentions(x, hidden::Vector{Symbol}) = false
    mentions(name::Symbol, hidden::Vector{Symbol}) = (name in names) && !(name in hidden)
    function mentions(e::Expr, hidden::Vector{Symbol})
        if e.head == :call
            any(a -> mentions(a, hidden), @view e.args[2:end]) && return true
            # A user function's body only sees its own parameters, not its caller's.
            f = e.args[1]
            if (f isa Symbol) && haskey(macros, f) && !(f in checked_macros)
                push!(checked_macros, f)
                m = macros[f]
                return mentions(m.body, m.params)
            end
            return false
        elseif e.head == :comparison
            # Operators are between the values.
            return any(i -> mentions(e.args[i], hidden), 1:2:length(e.args))
        else
            return any(a -> mentions(a, hidden), e.args)
        end
    end
    return mentions(p.expr, Symbol[ ])
end


###############
#  Compiling

mutable struct CalcCompiler
    env::CalcEnvironment
    input_names::Vector{Symbol}
    array_names::Vector{Symbol}
    # The environment which user functions' bodies are compiled in (plus their own parameters).
    # It's `env` itself, except within a user function's body, where `env` also has that function's parameters.
    macro_env::CalcEnvironment
end

"Gets the slot of the given input or array name, adding it if it's new"
function get_slot!(names::Vector{Symbol}, name::Symbol)::Int
    i = findfirst(==(name), names)
    if isnothing(i)
        push!(names, name)
        return length(names)
    else
        return i
    end
end

"
Checks and compiles the given Julia expression (or literal) into a calculation.
Throws a `CalcError` if it's not valid.
"
function calc_compile(expr, env::CalcEnvironment = CalcEnvironment())::CalcProgram
    return calc_compile(expr, env, env)
end
function calc_compile(expr, env::CalcEnvironment, macro_env::CalcEnvironment)::CalcProgram
    expr = normalize_expr(expr)
    compiler = CalcCompiler(env, Symbol[ ], Symbol[ ], macro_env)
    root = compile_node(compiler, expr)
    return CalcProgram(expr, root, compiler.input_names, compiler.array_names)
end

"
Copies the expression without any `LineNumberNode`s, and with single-expression blocks unwrapped,
  so that equivalent sources compare equal (e.g. `a ? b : c` and the `if` block Julia prints it as).
"
normalize_expr(x) = x
function normalize_expr(e::Expr)
    args = Any[ normalize_expr(a) for a in e.args if !isa(a, LineNumberNode) ]
    return ((e.head == :block) && (length(args) == 1)) ? args[1] : Expr(e.head, args...)
end

compile_node(::CalcCompiler, x::Union{Bool, Integer, AbstractFloat}) = CalcNodeConst(to_calc_value(x))
function compile_node(c::CalcCompiler, name::Symbol)
    if haskey(c.env.constants, name)
        return CalcNodeConst(c.env.constants[name])
    elseif name in c.env.inputs
        return CalcNodeInput(get_slot!(c.input_names, name))
    elseif haskey(c.env.functions, name) || haskey(c.env.array_functions, name)
        calc_error("`", name, "` is a function; call it like `", name, "(...)`")
    else
        calc_error("Unknown name `", name, "`")
    end
end
compile_node(::CalcCompiler, x) = calc_error("Unsupported value in a calculation: `", x, "` (a ", typeof(x), ")")

function compile_node(c::CalcCompiler, expr::Expr)
    if expr.head == :call
        return compile_call(c, expr)
    elseif expr.head == :comparison
        return compile_comparison(c, expr)
    elseif expr.head in (Symbol("&&"), Symbol("||"))
        @assert length(expr.args) == 2
        (a, b) = compile_node.(Ref(c), expr.args)
        node = (expr.head == Symbol("&&")) ? CalcNodeAnd(a, b) : CalcNodeOr(a, b)
        return try_fold(node, expr)
    elseif expr.head in (:if, :elseif)
        if length(expr.args) < 3
            calc_error("`if` needs an `else` branch, in `", expr, "`")
        end
        condition = compile_node(c, expr.args[1])
        if condition isa CalcNodeConst
            # Only compile the branch that will be taken.
            return compile_node(c, expr.args[expect_bool(condition.value, expr) ? 2 : 3])
        end
        return CalcNodeIf(condition, compile_node(c, expr.args[2]), compile_node(c, expr.args[3]))
    elseif expr.head == :block
        if length(expr.args) != 1
            calc_error("A block in a calculation must contain exactly one expression; got `", expr, "`")
        end
        return compile_node(c, expr.args[1])
    else
        calc_error("Unsupported syntax in a calculation (`", expr.head, "`): `", expr, "`")
    end
end

function compile_call(c::CalcCompiler, expr::Expr)
    (fn_name, raw_args...) = expr.args
    if !isa(fn_name, Symbol)
        calc_error("Unsupported function call in a calculation: `", expr, "`")
    elseif any(a -> Base.isexpr(a, (:parameters, :kw)), raw_args)
        calc_error("Keyword arguments aren't supported in a calculation: `", expr, "`")
    elseif any(a -> Base.isexpr(a, Symbol("...")), raw_args)
        calc_error("Splatting isn't supported in a calculation: `", expr, "`")
    end

    # Random interval:
    if fn_name == :(:)
        if length(raw_args) != 2
            calc_error("Random intervals must be `a:b` (no step size), got `", expr, "`")
        end
        (a, b) = compile_node.(Ref(c), raw_args)
        if (a isa CalcNodeConst) && (b isa CalcNodeConst)
            check_interval_types(a.value, b.value)
        end
        return CalcNodeInterval(a, b)
    end

    # Curved random interval:
    if fn_name == :curve
        if (length(raw_args) != 2) || !Base.isexpr(raw_args[1], :call) || (raw_args[1].args[1] != :(:))
            calc_error("Expected `curve(a:b, power)`, got `", expr, "`")
        end
        interval = compile_node(c, raw_args[1])::CalcNodeInterval
        power = compile_node(c, raw_args[2])
        if power isa CalcNodeConst
            check_curve_power(power.value)
        end
        return CalcNodeCurvedInterval(interval.a, interval.b, power)
    end

    # Rounding to Int, e.g. `floor(Int, x)`:
    if haskey(BUILTIN_ROUNDING_TO_INT, fn_name) && (length(raw_args) == 2) && (raw_args[1] == :Int) &&
       (c.env.functions[fn_name] === BUILTIN_FUNCTIONS[fn_name])
        return try_fold(CalcNodeCall(fn_name, BUILTIN_ROUNDING_TO_INT[fn_name], [ compile_node(c, raw_args[2]) ]),
                        expr)
    end

    # Array functions:
    if haskey(c.env.array_functions, fn_name)
        (array_name, query) = c.env.array_functions[fn_name]
        slot = get_slot!(c.array_names, array_name)
        if isnothing(query)
            if length(raw_args) != 1
                calc_error("`", fn_name, "(axis)` takes 1 argument; got `", expr, "`")
            end
            return CalcNodeArrayElement(fn_name, slot, compile_node(c, raw_args[1]))
        else
            if !isempty(raw_args)
                calc_error("`", fn_name, "()` takes no arguments; got `", expr, "`")
            end
            return CalcNodeArrayQuery(fn_name, slot, query)
        end
    end

    # User functions, if allowed, are anything else; they're looked up when the call runs.
    # The arguments belong to this calculation, so they're compiled now.
    if !haskey(c.env.functions, fn_name) && c.env.macros
        return CalcNodeMacroCall(fn_name, CalcNode[ compile_node(c, a) for a in raw_args ],
                                 collect(Any, raw_args), c.macro_env, Ref{Any}(nothing))
    end

    # Normal functions:
    if !haskey(c.env.functions, fn_name)
        calc_error("Unknown function `", fn_name, "`")
    end
    f = c.env.functions[fn_name]
    if !in(length(raw_args), f.n_args)
        calc_error("`", fn_name, "` takes ", format_arg_count(f.n_args), " arguments; got `", expr, "`")
    end
    return try_fold(CalcNodeCall(fn_name, f, CalcNode[ compile_node(c, a) for a in raw_args ]), expr)
end

function compile_comparison(c::CalcCompiler, expr::Expr)
    # Args alternate between values and operators: a < b <= c.
    args = CalcNode[ compile_node(c, a) for a in expr.args[1:2:end] ]
    op_names = Symbol[ ]
    ops = CalcFunction[ ]
    for op in expr.args[2:2:end]
        if !isa(op, Symbol) || !haskey(c.env.functions, op) || !in(2, c.env.functions[op].n_args)
            calc_error("Unsupported comparison `", op, "` in `", expr, "`")
        end
        push!(op_names, op)
        push!(ops, c.env.functions[op])
    end
    return try_fold(CalcNodeComparison(args, op_names, ops), expr)
end

children(n::CalcNodeCall) = n.args
children(n::CalcNodeComparison) = n.args
children(n::Union{CalcNodeAnd, CalcNodeOr}) = (n.a, n.b)

"
If all the node's children are constant, evaluates it now and returns a constant node.
If the evaluation fails, the node is kept as-is, so the error only happens if it's actually evaluated
  (e.g. in an `if` branch that's never taken).
"
function try_fold(node::CalcNode, expr)::CalcNode
    if all(n -> n isa CalcNodeConst, children(node))
        try
            return CalcNodeConst(eval_node(node, nothing))
        catch e
            (e isa CalcError) || rethrow()
        end
    end
    return node
end


################
#  Evaluating

# Deliberately not specialized on the RNG type, so that the precompiled evaluator works for any RNG.
struct CalcEvalState
    rng::AbstractRNG
    inputs::Vector{CalcValue}
    arrays::Vector{Vector{CalcValue}}

    # For user functions, which compile into their own programs that need their own inputs:
    macros::CalcMacros
    given::Any # Every input and array given to `calc_evaluate()`
    macro_depth::Int
end

"The default error handler for `calc_evaluate()`, which simply throws the error"
calc_throw(e::CalcError) = throw(e)

"
Evaluates the calculation, using the given RNG for random intervals.
Pass the values of the environment's inputs and arrays as keyword arguments;
  only the ones the calculation actually uses are required, and any others are ignored.

If something goes wrong, a `CalcError` (whose message names the calculation) is passed to `on_error`,
  whose return value becomes the result of the calculation.
By default it throws the error; you could instead log it and return a sane fallback value.

`macros` holds the user functions (see `CalcMacro`) which calls are looked up in,
  if the calculation was compiled in an environment that allows them.
"
function calc_evaluate(p::CalcProgram, rng::AbstractRNG = Random.default_rng(), on_error = calc_throw,
                       macros::CalcMacros = NO_CALC_MACROS
                       ; inputs...)::CalcValue
    if p.root isa CalcNodeConst
        return p.root.value
    end

    try
        return evaluate_program(p, rng, macros, values(inputs), 0)
    catch e
        (e isa CalcError) || rethrow()
        return to_calc_value(on_error(CalcError(string("In `", p.expr, "`: ", e.msg))),
                             "the error handler")
    end
end

"Evaluates a program, throwing any `CalcError` as-is"
function evaluate_program(p::CalcProgram, rng::AbstractRNG, macros::CalcMacros, given, macro_depth::Int)::CalcValue
    get_given(name) = haskey(given, name) ? given[name] : calc_error("Missing input `", name, "`")
    state = CalcEvalState(
        rng,
        CalcValue[ to_calc_value(get_given(name), "input `$name`") for name in p.input_names ],
        Vector{CalcValue}[ to_calc_array(get_given(name), name) for name in p.array_names ],
        macros, given, macro_depth
    )
    return eval_node(p.root, state)
end

function expect_bool(x::CalcValue, context)::Bool
    if !isa(x, Bool)
        calc_error("Expected true or false, got `", x, "`, in `", context, "`")
    end
    return x
end

eval_node(n::CalcNodeConst, s)::CalcValue = n.value
eval_node(n::CalcNodeInput, s::CalcEvalState)::CalcValue = s.inputs[n.slot]
function eval_node(n::CalcNodeArrayElement, s::CalcEvalState)::CalcValue
    array = s.arrays[n.slot]
    axis = eval_node(n.axis, s)
    if !isa(axis, Int)
        calc_error("`", n.name, "(axis)` needs an Int axis; got `", axis, "`")
    elseif !in(axis, 1:length(array))
        calc_error("`", n.name, "(", axis, ")` is out of range; the array has ", length(array), " elements")
    end
    return array[axis]
end
function eval_node(n::CalcNodeArrayQuery, s::CalcEvalState)::CalcValue
    array = s.arrays[n.slot]
    if isempty(array)
        calc_error("`", n.name, "()` can't be used on an empty array")
    end
    return to_calc_value(n.query(array), "function `$(n.name)`")
end
function eval_node(n::CalcNodeCall, s)::CalcValue
    # Avoid allocating for the common argument counts.
    args = n.args
    if length(args) == 1
        return call_calc_function(n.f, n.name, eval_node(args[1], s))
    elseif length(args) == 2
        return call_calc_function(n.f, n.name, eval_node(args[1], s), eval_node(args[2], s))
    elseif length(args) == 3
        return call_calc_function(n.f, n.name, eval_node(args[1], s), eval_node(args[2], s),
                                               eval_node(args[3], s))
    else
        return call_calc_function(n.f, n.name, (eval_node(a, s) for a in args)...)
    end
end
function eval_node(n::CalcNodeComparison, s)::CalcValue
    a = eval_node(n.args[1], s)
    for i in 1:length(n.ops)
        b = eval_node(n.args[i + 1], s)
        result = call_calc_function(n.ops[i], n.op_names[i], a, b)
        if !expect_bool(result, n.op_names[i])
            return false
        end
        a = b
    end
    return true
end
eval_node(n::CalcNodeAnd, s)::CalcValue = expect_bool(eval_node(n.a, s), "&&") &&
                                          expect_bool(eval_node(n.b, s), "&&")
eval_node(n::CalcNodeOr, s)::CalcValue = expect_bool(eval_node(n.a, s), "||") ||
                                         expect_bool(eval_node(n.b, s), "||")
eval_node(n::CalcNodeIf, s)::CalcValue = if expect_bool(eval_node(n.condition, s), "if")
    eval_node(n.if_true, s)
else
    eval_node(n.if_false, s)
end
eval_node(n::CalcNodeInterval, s::CalcEvalState)::CalcValue =
    sample_interval(eval_node(n.a, s), eval_node(n.b, s), s.rng)

"Numbers which aren't `Bool`; an interval mixing `Int` and `Float64` is a `Float64` interval"
const CalcNumber = Union{Int, Float64}

check_interval_types(a::CalcNumber, b::CalcNumber) = nothing
check_interval_types(a::Bool, b::Bool) = nothing
check_interval_types(a, b) = calc_error(
    "A random interval can't mix `Bool` with numbers; ",
    "got `", a, "` (", typeof(a), ") and `", b, "` (", typeof(b), ")"
)

"Picks a random value between `a` and `b` inclusive, in either order"
sample_interval(a::Int, b::Int, rng::AbstractRNG)::Int = rand(rng, min(a, b):max(a, b))
function sample_interval(a::Float64, b::Float64, rng::AbstractRNG)::Float64
    (lo, hi) = minmax(a, b)
    return lo + ((hi - lo) * rand(rng, Float64))
end
sample_interval(a::CalcNumber, b::CalcNumber, rng::AbstractRNG)::Float64 =
    sample_interval(convert(Float64, a), convert(Float64, b), rng)
sample_interval(a::Bool, b::Bool, rng::AbstractRNG)::Bool = (a == b) ? a : rand(rng, Bool)
function sample_interval(a, b, rng)
    check_interval_types(a, b)
    error("Unhandled: ", typeof(a))
end

function eval_node(n::CalcNodeCurvedInterval, s::CalcEvalState)::CalcValue
    a = eval_node(n.a, s)
    b = eval_node(n.b, s)
    check_interval_types(a, b)
    power = check_curve_power(eval_node(n.power, s))
    return sample_curved(a, b, det_pow(rand(s.rng, Float64), power))
end

function eval_node(n::CalcNodeMacroCall, s::CalcEvalState)::CalcValue
    m = get(s.macros, n.name, nothing)
    if isnothing(m)
        calc_error("Unknown function `", n.name, "`")
    elseif length(m.params) != length(n.args)
        calc_error("User function `", n.name, "` takes ", length(m.params), " arguments; got ", length(n.args))
    elseif s.macro_depth >= MAX_MACRO_DEPTH
        calc_error("User function calls are nested more than ", MAX_MACRO_DEPTH, " deep, ",
                   "at `", n.name, "`; is it infinitely recursive?")
    end

    # Evaluate the arguments once each, before the body sees them.
    arg_values = ntuple(i -> eval_node(n.args[i], s), length(n.args))

    try
        # Compile the body if the macro is new or was replaced since last time.
        (cached_macro, program) = something(n.cache[], (nothing, nothing))
        if cached_macro !== m
            program = calc_compile(m.body, calc_macro_body_env(n.env, m.params), n.env)
            n.cache[] = (m, program)
        end
        # The body sees the caller's inputs, with the parameters on top.
        given = merge(s.given, NamedTuple{Tuple(m.params)}(arg_values))
        return evaluate_program(program, s.rng, s.macros, given, s.macro_depth + 1)
    catch e
        (e isa CalcError) || rethrow()
        calc_error("In user function `", n.name, "(", join(n.arg_exprs, ", "), ")`: ", e.msg)
    end
end

function check_curve_power(p::CalcValue)::Float64
    if (p isa Bool) || !(p > 0)
        calc_error("The power of a curved interval must be a number above 0; got `", p, "`")
    end
    return convert(Float64, p)
end

"Maps `t` (from 0 to 1) onto the interval, going from `a` toward `b`"
sample_curved(a::Float64, b::Float64, t::Float64)::Float64 = lerp(a, b, t)
sample_curved(a::CalcNumber, b::CalcNumber, t::Float64)::Float64 =
    sample_curved(convert(Float64, a), convert(Float64, b), t)
function sample_curved(a::Int, b::Int, t::Float64)::Int
    # Split into equal bins, one per value, starting from `a`.
    n_values = abs(b - a) + 1
    bin = min(floor(Int, t * n_values), n_values - 1)
    return a + (bin * sign(b - a))
end
sample_curved(a::Bool, b::Bool, t::Float64)::Bool = Bool(sample_curved(Int(a), Int(b), t))


#################
#  Precompiling

# Exercise every kind of node, so the first real calculation doesn't wait on the JIT.
@compile_workload begin
    env = CalcEnvironment(inputs = (:x, ), arrays = (:a, ), constants = (k = 2, ),
                          functions = (twice = CalcFunction(v -> 2v, 1), ))
    for src in (
        :( clamp(x * k + a(1), 0, 10) / twice(a_max()) ),
        :( 0 < x <= a_sum() && !(x == a_min_axis()) || a_median() > a_average() ),
        :( x > 2 ? (1:6) + (0.5:1.5) + curve(9:1, x) : floor(Int, sqrt(x)) ),
        :( curve(false:true, 2) ? curve(0.0:float(x), 0.5) : x ),
        :( if false:true
               lerp(0.0, a_prod(), inv_lerp(1, a_len(), x))
           else
               max(a_min(), x, a_max_axis(), exp(-x) ^ 1.5, log(2, x + 1), log10(x) + exp2(-x))
           end )
    )
        program = calc_compile(src, env)
        calc_evaluate(program, Random.Xoshiro(1); x = 3, a = (1, 2, 3))
        calc_evaluate(program, Random.Xoshiro(1), e -> 0; x = 3.5, a = [ 1.5, 2.5 ])
        string(program)
    end

    macro_env = CalcEnvironment(inputs = (:x, ), arrays = (:a, ), macros = true)
    macros = Dict(calc_macro_parse(:( half(v) = v / 2 )),
                  calc_macro_parse("k() = half(a_sum()) + 1"))
    program = calc_compile(:( half(x) * k() ), macro_env)
    calc_evaluate(program, Random.Xoshiro(1), calc_throw, macros; x = 3, a = (1, 2, 3))
    string(calc_macro_expr(:half, macros[:half]))
end

end # module
