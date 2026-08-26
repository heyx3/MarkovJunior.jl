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

Read-write data at runtime, for one specific run of an algorithm,
  is kept in the `data_store` dictionary:

````julia
# (inside the implementation of a custom Op or Bias)
algo_state.data_store[:my_custom_n_applications] += 1
depth_i::Int = algo_state.data_store[:my_custom_n_applications]
````

> *Keep in mind `data_store` is a type-unstable `Dict{Symbol, Any}`,*
> *so make sure the data you pull from it is annotated with a type!*
> *You may even want to wrap the value in a `RefValue` so it can be mutated without boxing.*

## Allocators

Users can provide their own "allocators" to control the many arrays and sets that are used
  in each run of an algorithm.
This allocator is provided in the algorithm-state object,
  and you should use it when creating new collections.
The allocator itself is type-unstable, so you need to provide type annotations for its output!

````julia
# Our custom Op needs to grab a new Int32 array at half the grid resolution.
my_alloc::Array{Int32, ndims(grid)} = markov_allocator_acquire_array(
    algo_state.allocator,
    size(grid) .÷ 2,   Int32
)
...
# Our Op is finished with the array.
markov_allocator_release_array(algo_state.allocator, my_alloc)
````

For a resizable `Vector` that starts empty, pass a 0-tuple for the array size.

You can also allocate `Set`s, and 'ordered sets' from the package *OrderedCollections.jl*,
  in a similar way but without providing a size.
Each new allocated set is expected to already be empty.

> *OrderedSets are important because all algorithm logic must be deterministic!*
> *Never iterate through normal `Set` and `Dict` objects.*

## Logging

There is a bult-in system of logging macros for algorithm logic.
When implementing your own Ops or Biases, you should throw data into it!

It is controlled by a "compile-time" flag, which in Julia can actually be toggled at will --
  albeit at the cost of significant JIT overhead right after toggling.
It defaults to off, and you can turn the logging on by redefining the function `MarkovJunior.log_logic() = true`.

The following logger macros are available (but un-exported):

* `@logic_log(...)` and `@logic_logln(...)` will print the given arguments, with a tab inserted after every newline.
* `@logic_tab_in()` will increase the tab level.
* `@logic_tab_out()` will decrease the tab level.

## Parsing

Extensions to the DSL must be able to parse themselves from a user string,
  and convert themselves back into string form.
You are not expected to preserve cosmetic details like comments or whitespace.

Parsing is done in terms of a standard Julia syntax tree ("AST"),
  meaning the parsed code must be made of valid Julia syntax structures.
It does *not* necessarily have to be made of valid Julia code.
For example `4 = 5` is easily parsed, so it's allowed within our DSL,
  but would of course throw an error if directly compiled and executed.
Unless your custom Op syntax is trivial, use *MacroTools.jl* to simplify searching and transforming the AST.

While parsing your new extensions you will have access to a "stack trace",
  in the form of two helper functions.
This trace can report exactly where errors occur in the DSL.
To use it, you only need two parameters passed into your parsing fuction, `loc` and `inputs`:

* `raise_parse_error(loc, inputs, msg...)` will raise an error at the current parsing stack trace.
  Always use this instead of `error` (unless you encounter a genuine internal error).
* `with_parser_stacktrace(inputs, msg_str) do ... end` will push a new element to the stack trace,
  execute your block of code, then pop back off the stack trace.
  It is recommended to wrap your code in this any time you go deeper into your DSL expressions.

The stack trace is stored in `inputs.op_stack_trace::Stack{Any}`, in case you want to manipulate it yourself.

## Custom Ops

If the likes of `@fill`, `@rewrite`, and `@sequence` aren't enough for you,
  you can define your own grid operations with their own DSL syntax.

First create a struct inheriting from `AbstractMarkovOp`,
  with whatever data is needed to define your Op.
The data in this struct should closely match the features of your planned DSL syntax.

Now simply implement `markov_algo_run` for your Op:

````julia
const MJ = MarkovJunior
function MJ.markov_algo_run(
    op::MyCustomOp,
    algo::MarkovAlgorithm,
    # Contains the grid, RNG, allocator, data_store, etc.
    algo_state::MJ.AlgoState,

    # Some number of pre-existing biases may be active when your Op executes:
    inherited_biases::Tuple{Vararg{MJ.AbstractMarkovBias}},
    inherited_bias_states::Tuple{Vararg{Any}},

    # Often you will want more info at compile-time, such as the exact grid type.
    # The best way to do that is with hidden extra arguments:
    grid::TGrid = algo_state.grid,
    alloc::TAllocator = algo_state.allocator
)::Tuple{
    # Your function returns `(changed_anything, new_inherited_bias_states)`
    Bool, typeof(inherited_bias_states)
} where {
    NGridDims, TGrid::MJ.CellGid{NGridDims},
    TAllocator<:MarkovJunior.AbstractMarkovAllocator
}
    VI = Vec{Int32, NGridDims} # typedef for grid index

    # Some Ops want to reallocate the grid with a new size or dimensionality.
    # The proper way to do this is to call our helper function:
    grid = MJ.markov_algo_new_grid(algo_state, size(grid) .÷ 2) do new_grid, old_grid
        # (fill in 'new_grid', referencing 'old_grid')
        return nothing
    end

    # Inherited biases need to be notified as you make changes to the grid.
    # You must pass the specific rectangular area you modified, and the grid's previous values there.
    # If your op is very broad, or you want to be lazy,
    #   just wait till the end and say that the whole grid changed:
    return MJ.markov_allocator_with_array(alloc, size(grid), eltype(grid)) do old_grid
        old_grid .= grid

        # (do stuff to the grid)

        # Every so often, emit tick events of various priorities.
        # Tick calls are blocking, until the user acknowledges them
        #    and decides the algorithm should resume.
        # In other words, this algorithm is a coroutine!
        #
        # The lowest priorities (1 and 2) are for extremely minor events,
        #    like changing individual pixels, and usually get compiled out for performance.
        #
        # The min tick priority that doesn't get compiled out is
        #    `MJ.STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY`, currently 3.
        # The standard priority for the *completion* of a simple Op like @rewrite is
        #    `MJ.STANDARD_END_OF_OP_TICK_PRIORITY`, currently 4.
        MJ.markov_algo_tick(algo_state, MJ.STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY)

        # You can also tick with a Symbol, which represents a "tagged event"
        #    that gets sent to users.
        # This is what the `@event` Op does.
        MJ.markov_algo_tick(algo_state, :my_custom_event)
        # One use of tagged events is to let users inject their own changes to the grid,
        #    then your Op can react to those changes.
        # They could also feed more complex data to you through `algo_state.data_store`.

        # Here's some sample Op logic:
        my_area = ...
        for pixel::VI in my_area
            grid[pixel] = 0
            MJ.markov_algo_tick(algo_state, 1)
        end
        MJ.markov_algo_tick(algo_state, MJ.STANDARD_END_OF_OP_TICK_PRIORITY)

        # Lastly, we need to report any changes and update the Biases to know about our changes.
        return if isempty(my_area)
            (false, inherited_bias_states)
        else
            (true, MJ.markov_bias_update.(
                inherited_biases,  inherited_bias_states,
                algo, algo_state,
                # Our update changed the entire grid, more or less
                BoxI{NGridDims}(min=one(VI), size=vsize(grid)),
                old_grid
            ))
        end
    end
end
````

> *Practice good exception hygiene: always free your allocated collections,*
> *and make sure that you don't stop `MarkovJunior.ErrorCancelAlgo` from propagating upward!*

Lastly, implement parsing and string-ification.

````julia
using MacroTools # Very helpful to parse Julia expressions, using @capture

# Suppose your Op's syntax is `@my_op (a, b) c...`
MJ.dsl_string(op::MyCustomOp) = let a = op.a,
                                    b = op.b,
                                    cs = join(("($c)" for c in op.c), " ")
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

    cs = [ ]
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

Now you can use this new Op in your algorithms!

## Custom Biases

Biases influence the `@rewrite` Op to prefer some rewrites over others.
They can also be added to a `@sequence`,
  causing them to be inherited by every Op within that sequence.
New custom Ops may choose to use biases as well.

The workflow for new Biases is more complicated than new Ops:

1. "Initialize" the Bias and return the initial state: `markov_bias_initialize`
2. "Update" the Bias after a grid modification, and return its new state: `markov_bias_update`
3. Get the bias value for a potential change to the grid: `markov_bias_calculate`
  The output must be a non-negative number, or `nothing` if the change is strictly forbidden.
4. "Cleanup" the Bias's allocations: `markov_bias_cleanup`
5. Convert the Bias to and from a DSL string: `dsl_string`, `parse_markovjunior_bias`
  The DSL syntax for a bias is a function call (`my_bias(a, b/c)`),
  compared to a macro for Ops (`@my_op a b/c`).
1. Add constraints to where and how the bias is used (e.g. only one may exist at a time):
  `markov_bias_validate`.
  See its doc-string for more info on how it works.

## Custom Priorities

Priorities decide which rule a `@rewrite` Op focuses on first.
They are similar to biases, but live within a single Op
  and work per-rule rather than per-application of a rule.
The DSL syntax is `PRIORITIZE(name, args...)`, for example `PRIORITIZE(earliest)`.

**NOTE: the interface for priorities will change in the future, so be careful when upgrading MarkovJunior.jl**

1. Define a struct inheriting from `AbstractMarkovRewriteProperty`.
2. Define its conversion to/from string: `dsl_string` and `parse_markovjunior_rewrite_priority`
   * For more info on how to parse things, see the parsing section for custom Ops above.
3. Add the priority's logic by implementing `pick_rule_using_rewrite_priority`.

See the bulit-in priorities for working examples.