# Extending MarkovJunior.jl

There are multiple ways to extend this package.
To start, understand the mechanisms for providing custom data to the algorithm.

## Custom Data

### Pragmas

To embed data in the `@markovjunior` language itself,
  add `@pragma` statements to the top of the main block.

````julia
mj = @markovjunior begin
    @pragma A b c d
    @pragma B
    @pragma(A, 1, 2, 3) # Julia has two syntaxes for macros; this is the other one

    ... # Normal algorithm stuff goes here
end

# Pragmas are normally accessed from a lookup of Name to Args.
# Each list of args is stored in chronological order.
@assert(mj.pragmas_map[:A] == [ [ :b, :c, :d ], [ 1, 2, 3 ] ])
# You can also look up the chronological order of EVERY pragma statement.
@assert(mj.pragmas_chronological == [ :A => 1, :B =>  1, :A => 2 ])
````

### Add-ons

To add parameters to a specific parsed `@markovjunior` instance, write into its `add_ons` dictionary.
This allows you to set data per-instance, even when parsing from the same algorithm string.

````julia
mj = @markovjunior ... # Some algorithm
mj.add_ons[:my_custom_op_param] = 4.5  # Only this instance has this param
````

> *You should not modify `add_ons` while an instance of the algorithm is running!*
> *For data that changes during a run, see the next section.*

### Data Store

When extending the language, you can store and modify data at runtime.
Use the `data_store` dictionary, in the algorithm state that's passed to you:

````julia
# (inside the implementation of a custom Op or Bias)
state.data_store[:my_custom_n_applications] += 1
depth_i::Int = state.data_store[:my_custom_n_applications]
````

**NOTE**: because `data_store` is a type-unstable `Dict{Symbol, Any}`,
  I highly recommend storing your custom data in a mutable form (like `Ref{T}`)
  which you can pull into a type-stable local variable.

## Allocators

Users can provide their own "allocators" to control the many arrays and sets that are used
  in each run of an algorithm.
This allocator is stored in the algorithm state, available to your custom extensions.

````julia
# Our custom Op needs to grab a new Int32 array at half the grid resolution.
# The allocator is type-unstable, so it's important to specify the type of the variable.
my_alloc::Array{Int32, ndims(grid)} = markov_allocator_acquire_array(
    state.allocator,
    size(grid) .÷ 2, Int32
)
# If you're grabbing several allocations at once, dispatch on the allocator type to make this easier:
function with_alloc(alloc::TAlloc) where {TAlloc}
    # Now allocator calls are type-stable!
    ...
end
with_alloc(state.allocator)
...
# Our Op is finished with the array.
markov_allocator_release_array(state.allocator, my_alloc)
````

Most of the time you should use our scoped wrapper:

````julia
markov_allocator_with_array(state.allocator, size(grid) .÷ 2, Int32) do small_grid
    ... # Use the grid
end
````

You can also allocate sets and 'ordered sets' (from the package *OrderedCollections.jl*),
  in a similar way but without providing a size.
You can expect the sets from the allocator to always come empty.

> *OrderedSets are important because all algorithm logic must be deterministic!*
> *Don't ever iterate through normal `Set` and `Dict` objects.*

## Logging

There is a built-in system of macros for logging algorithm logic.
When implementng your own Ops or Biases, you should throw data into it!

It is controlled by a "compile-time" flag, which in Julia can actually be toggled at will --
  albeit at the cost of some JIT compilation time.
It defaults to off, and you can turn the logging on by redefining the function `MarkovJunior.log_logic() = true`.

The following logger macros are available (but un-exported):

* `@logic_log(...)` and `@logic_logln(...)` will print the given arguments, with a tab inserted after every newline.
* `@logic_tab_in()` will increase the tab level.
* `@logic_tab_out()` will decrease the tab level.

## Parsing

Most extensions to the DSL must be able to convert themselves into their string form.
Obviously they also must be able to parse themselves *from* a string,
  more specifically a Julia syntax tree ("AST") that comes from having Julia parse the string.
If your custom Op syntax is non-trivial, use *MacroTools.jl* to simplify searching and transforming the AST.

Within parsing functions you will also have access to a MarkovJunior parser "stack trace",
  in the form of two helper functions.
These functions need some of the parameters passed into your parsing function, namely `loc` and `inputs`.

* `raise_parse_error(loc, inputs, msg...)` will raise an error that includes the current parsing stack trace.
* `with_parser_stacktrace(inputs, msg_str) do ... end` will execute a block of code
with your custom message pushed onto the parsing stack trace.

The stack trace is stored in `inputs.op_stack_trace::Stack{Any}`, in case you want to manipulate it yourself.

## Custom Ops

If the likes of `@fill`, `@rewrite`, and `@sequence` aren't enough for you,
  you can define your own grid operations.

First create a struct inheriting from `AbstractMarkovOp` defining what your Op will do.

**TODO: Finish**

If you want your Op to reallocate the grid, call `markov_algo_new_grid(state, size_tuple)::CellGrid`.
The algorithm state's `grid` field will point to the new array after this call.

Lastly, implement parsing and string-ification of your Op.

````julia
using MacroTools # Very helpful to parse Julia expressions, using @capture

# Suppose your Op's syntax is `@my_op (a, b) c...`
MJ.dsl_string(op::MyCustomOp) = let a = op.a,
                                    b = op.b,
                                    cs = join(string.(op.c), " ")
    "@my_op ($a, $b) $cs"
end
function MJ.parse_markovjunior_op(::Val{Symbol("@my_op")}, inputs::MJ.MacroParserInputs,
                                  loc::LineNumberNode,
                                  expr_args,
                                  full_line::Expr)
    if length(expr_args) < 1
        MJ.raise_parse_error(loc, inputs, "Expected at least one argument in the form `(a, b)`")
    elseif !@capture(expr_args[1], (a_, b_)) || !isa(a, Integer) || !isa(b, Integer)
        MJ.raise_parse_error(loc, inputs,
                             "First argument should be a tuple of two integers; got ", expr_args[1])
    end

    cs = Int[ ]
    for (c,i) in enumerate(expr_args[2:end])
        MJ.with_parser_stacktrace(inputs, "c[$i]: $c") do
            if !isa(c, Integer)
                MJ.raise_parse_error(loc, inputs, "Expected integer!")
            else
                push!(cs, convert(Int, c))
            end
        end
    end

    return MyCustomOp(convert(Int, a), convert(Int, b), cs)
end
````

## Custom Biases

Biases influence the `@rewrite` Op to prefer some rewrites over others.
They can also be added to a `@sequence`, causing them to be inherited by every Op within that sequence.
The workflow for Biases is almost identical to Ops:

**TODO: Rewrite, and then update the below ways in which it's different from Ops

* The "update" happens once per rewrite operation, allowing it to react to any changes to the grid
* Between each update, the Bias has a "calculate" function which is called many times --
to decide the bias towards each potential rewrite move: `markov_bias_calculate`
  * The output of calculate must be a non-negative number, or `nothing` if the move is strictly illegal.
* The DSL syntax for a bias is a function call (`my_bias(a, b/c)`) rather than a macro (`@my_op a b/c`)
* After being parsed the Bias can enforce constraints,
  for example you can forbid more than one of itself to exist on a single Op.
  This is done with `check_markovjunior_biases`.

## Custom Priorities

Priorities decide which rule a `@rewrite` Op focuses on first.
They are similar to biases, but operate per-rule rather than per-application of a rule at a specific place.
The DSL syntax is `PRIORITIZE(name, args...)`, for example `PRIORITIZE(earliest)`.

1. Define a struct inheriting from `AbstractMarkovRewriteProperty`.
2. Define its conversion to/from string: `dsl_string` and `parse_markovjunior_rewrite_priority`
   * For more info on how to parse things, see the parsing section for custom Ops above.
3. Add the priority's logic by implementing `pick_rule_using_rewrite_priority`.