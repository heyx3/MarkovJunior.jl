# `@markovjunior` DSL syntax

Below is a comprehensive reference for the syntax to define a MarkovJunior algorithm.
If you're starting out, I recommend reading the example scenes and only coming here to answer specific questions.

Our new Julia macro `@markovjunior [dims] [clear_value] begin ... end`
  evaluates into an instance of `MarkovAlgorithm`.
This represents a sequence of operations that generate a grid of colored pixels,
  in any dimension (most commonly 2 for images and 3 for voxel scenes).
You could also think of it as defining an *animation*,
  which transforms a blank grid into an image/scene.

You can turn a parsed algorithm back into a DSL string with `dsl_string(algo)`.
However the result looks much worse than the original -- whitespace, comments and other niceties are lost!

## Notes about Julia

The DSL parser is written on top of Julia's own built-in ability to parse Julia code.
DSL statements do not have to be valid compilable Julia code,
  but they DO have to be made of parseable Julia syntax structures.
For example `4 = 5` is absurd code but totally valid to the parser,
  so it could in theory show up in our DSL.
Some operations, especially `@rewrite`, stretch this to the limit with pretty tortured syntax,
  but it's all parseable Julia at the end of the day.

For convenience, Julia allows macros (code statements with the `@` symbol) to be written two ways:

* Without parentheses and commas, to make statements shorter: `@my_macro a b c`
* With parentheses and commas, like a function call, to make multiline statements easier: `@my_macro(a, b, c)`.

If using the former we have a few inner expressions you *must* wrap in parentheses,
  which we will mention as they come up.

> ***Important note**: Julia uses 1-based indices, and so does our library! For example, the first axis is 1 and not 0.*

## Grid pixels 

Each pixel is a 1-byte uint storing a bitmapped color.
Types/names of colors are specified directly in the codebase.

Each color also has a single-char name which is used in this DSL, for example `R` stands for Red.

## Parameters

By default the grid starts with all pixels black, but you can change this initial color.
For example `@markovjunior 'M' begin ... end` starts with a Magenta fill.

> *If you want to customize the starting state, such as placing one red pixel in the center of the grid,*
> *use the [`@fill` operation](#fill).*

By default the grid can have any number of dimensions, but you can fix this with
  e.g. `@markovjunior 2 begin ... end`.
There are also some expressions inside the macro which imply a minimum number of dimensions,
  e.g. `@fill` with a 2D box will require at least a 2D grid.
If there are any dimensional mismatches between the `@markovjunior` and the grid, an error is thrown.

## Main Sequence

Inside the `begin ... end` block of `@markovjunior` should be a chronological sequence of operations,
  optionally preceded by `@pragma` notes.
For example:

````julia
@markovjunior begin
    @pragma viz 2     # Special notes to the interpreter, no effect on the core algorithm
    @sequence (length/12) begin    # Do this a few times, proportional to grid length
        @rewrite 1  b => R         # Place a single red pixel somewhere in the image
        @rewrite Rbb => RwR        # Randomly connect new pathways to the red pixel
        @rewrite w => R            # Finish the pathways to make a perfect maze
        @rewrite area/100  b => R  # Randomly remove some maze walls to form rooms, proportional to grid volume
        @rewrite Rb => RG          # Add walls around this area before starting the next one.
    end
    # Now open some trivial connections between walls, if possible.
    @rewrite (area/100) begin
        wRGGRw   => wwwwww
        wRG_GRw  => wwwwwww # Stomp over whatever was between the two walls, to make things interesting
        wRG__GRw => wwwwwwww # Slightly larger stomp
    end

    # Emit a "tagged event", which pauses the algorithm for users to view its current state.
    @event BeforeFinalSteps
    # If users need to *modify* the grid during a tagged event, prepend the name with '!'.
    # There's no special checking around this, so it's up to you to be vigilant!
    @event !BeforeFinalSteps

    # Emit another tagged event which gives users 
    # Finally, place a "player start" square.
    @fill B pixel(min=0, size=4)
    @fill T pixel(min=0, size=1)
end
````

## `@pragma`

This is [one of the standard ways](add_ons.md) to add external data to an algorithm or set algorithm-wide settings
**Note** that they are always specified at the top of the algorithm.
If you want something more like a live feature that can change during algorithm runs,
  implement a [custom Op](add_ons.md#custom-ops) that writes to [the live data-store](add_ons.md#data-store).

Our built-in GUI tool checks pragmas for the 3D materials of different cell types:

````julia
# Red should be desaturated stained-glass.
# The data format for glass is "opacity  roughness  albedo".
@pragma GuiMaterial R glass 0.5  0.2  (1.0, 0.5, 0.5)
# Blue should be empty space.
@pragma GuiMaterial B empty_space

# Green should be shiny, saturated plastic.
# The data format for non-metallic opaque surfaces is "roughness  albedo".
@pragma GuiMaterial G dielectric  0.2 (0.05, 1.0, 0.02)
#  (metal has the same format as dielectric)

# Yellow should be a bright light source.
# The data format for light sources is just the emissive HDR color.
@pragma GuiMaterial Y light_source  (5.0, 5.0, 1.0)
````

As another example, the `@fill` Op looks for `@pragma fast_fills`
  and, if it finds it, will always finish itself in a single tick.

## `@event`

This is a no-op which inserts a "tagged event", sort of like a named tick.
It signals to users that the algorithm has reached a certain point,
   and those users may wish to read the grid state before the algorithm continues to modify it.

````julia
@event BeforeSpawnPlacement
````

When done this way, users should *not* modify the grid! Only read from it.
To support modification, prepend the event name with `!`:

````julia
@event !BeforeSpawnPlacement
````

This will give secondary state, like all active [Biases](#bias-and-weights),
  a chance to update themselves after your modifications.
It can have a big performance impact under really active Biases, so use it sparingly.

## `@rewrite`

The primary operation in this algorithm is `@rewrite`, which executes one or more rewrite rules
  until no more rule applications are possible or some threshold is reached.
The rules are applied in a uniform-random order,
  though you can [use weights and "bias" terms](#bias-and-weights) to change that.

The official syntax for this operation is `@rewrite [threshold] rules [bias]`;
  the rules can be a single rule or a `begin ... end` block of them;
  the syntax for one such rule is `source => dest [modifiers...]`.
We'll go into detail on everything, but here is a quick cheat sheet of **all** the features:

````julia
@rewrite  (length/2 : area/50) #= <-- a Threshold range for how many times to run =#     begin
    # Rules go here.
    # If you only have one, you don't need to wrap it in this `begin end` block!

    # You may set a Priority before providing any rules.
    # These have higher precedence than anything else, including biases and weights.
    PRIORITIZE(earliest)

    R => b        # Basic rule
    R => G   * 2  # Twice as likely as other rules
    R => G  %0.75 # Forbidden in exactly 25% of grid pixels, randomly chosen when the @rewrite starts.
    R => G  %(0.5:0.9) # Forbidden in anywhere from 10% to 50% of grid pixels, randomly chosen when the @rewrite starts.

    # Underscore means "wildcard" on source pixels and "don't change" on destination pixels.
    R_B => __G

    YY[Gb] => RYY     # Third source pixel may be Green or black
    YY[Gb] => YY[bG]  # Third pixel swaps between Green and black
    YY[Gb] => YY{Gb}  # Third pixel may or may not swap

    Y__    => Y[3][2]   # Second and third pixel swap places

    # This rule only matches in the +X and +/-Y directions.
    # It also restricts the algorithm to 2D grids and higher, since it mentions a second axis.
    RGB => UMb  \[+x, Y]
    # You could also write those symmetries using numbers.
    # With letters you only have XYZW, so beyond 4D you need to use these,
    #    however to pick a single direction you also have to put the number in parentheses.
    # Add 'W...' or '4...' to the end to let it lay along all subsequent dimensions starting at 4D,
    #    in both directions.
    RGB => YMb  \[+(1), 2, 4...]

    # Order of modifiers is important -- write the symmetry *after* the weighting!
    # You can also divide the weight instead of multiply.
    RGB => YMb  /2  \[+x]

    # Get ready: here's a wacky rule that uses all the features at once.
    R_[Bb]w => [2]_[bB]{wbR}  %(0.4:0.6)  *4   \[-(1), z...]

    # One last thing: multidimensional rewrite rules (with multidimensional symmetries)!
    # They can do everything the above rules can do, but I'm keeping it simple in this example.
    [
        R G B
        w b g ;;;
        # ^^ One Z-slice, 'RGB' is along the first axis and 'Rw' the second
        # vv Another Z-slice, 'RS' is along the third axis
        S T U
        O w w ;;;;
        # Now start a new 3D slice along the fourth dimension!
        # 'RR' is along the fourth axis.
        R R R
        R R R ;;;
        R R R
        R R R
    ] => [ # Above was the source, below is the destination:
        R R R
        R R R ;;; # New Z-slice
        G G G
        G G G ;;;; # New 3D slice
        B B B
        B B B ;;;
        b b b
        b b b
    ] \[ # Now the symmetry modifier: Allow the block to only flip along the Z axis and swap the X/Y axes.
        (x, y)[ (+x, +y), (+y, +x) ]
        # Z is the only choice left for the block's Z,
        #   and not specifing anything means it can flip either way along that Z axis.
        # However we can add a "chirality" constraint between it and the Y axis,
        #   so that their handedness is preserved.
        {y, z}
    ]
end begin
    # Biases go here (see below).
    # If you only have one, you don't need to put it in a `begin end` block!
    field(G <- R, soft)
    field(G -> Y, recompute)
end
````

### Rewrite rules

The format of a single rewrite rule is `source => dest [modifiers]`,
  where both `source` and `dest` are an ordered list of color characters.
These lists represent a strip of pixels on the grid;
  the algorithm finds instances of `source` and potentially replaces them with `dest`.

The different modifiers are explained in detail below, but here is a quick reference in the same order they must appear in:
  * The *mask* `%X` randomly forbids a specific amount of the grid from matching the first pixel of this rule,
for example `%0.75` will randomly pick about one-fourth of the grid's pixels to not have the rule start there.
You can also provide a randomized range in parentheses, like `%(0.5:0.9)`.
Importantly, **all** rules in a single `@rewrite` op will use the same grid mask,
for example a rule with `%0.2` covers a strict subset of the area affected by a rule with `%0.5`.
  * The *weight* `*X` and `/X` change the chance of this rule being applied, relative to the others.
For example `*2` makes the rule twice as likely to be chosen.
  * The *symmetry* `\[ ... ]` lists the axes and directions this strip can run along:
`\[ x,  -y ]` can run along -X, +X, and -Y.
Add an ellipsis to allow every symmetry past a certain dimension, e.g. `W...` 
  allows the rule to face down any axis from W up to infinity.

Each color list is a string of color chars, e.g. `RGB` means "Red then Green then Blue".
Source and destination strings must describe the same number of pixels,
  so `RGB => wbgY` would throw an error.

You may provide a set of colors for a single source pixel,
  e.g. `RG[Bb]` means "Red then Green then (Blue or black)".
In this case you may make a corresponding list in the destination pixel, e.g.
  `RG[Bb] => RG[Yw]` means to replace the strip with "Red then Green then (Yellow if we were Blue; white if we were black)".
Destination color sets must have as many elements as their corresponding source.

You can use underscores, as a wildcard for source pixels and as "leave unchanged" for dest pixels.
For example, the rule `R_[Bb] => BY_` means to find a strip of "Red then anything then (Blue or black)",
  and replace with "Blue then Yellow then (whatever was there already)".

Destination colors can take randomly from a set, by enclosing their options with braces:
  `RGB => RG{wb}` turns Blue into either white or black.

Finally, destination colors can take specific colors from their source by index:
  `R_R => RR[2]` moves any color sandwiched between Reds to their edges.

### Multiple rules

You can provide a group of rewrite rules in a `begin ... end` block, for example:

````julia
@rewrite begin
    # All Red pixels become either Green or white.
    R => G
    R => w
    # All Yellow-Black pairs turn the black into white or Green.
    Yb => _[wG]
end
````

By default they all have an equal chance of being selected
  (weighted by count, meaning if there's twice as many candidates for the second rule then it's twice as likely to be chosen).
To change this, see [Priority](#priorities) and [Bias and weights](#bias-and-weights) below.

### Priority

In a block with multiple rules, you can change which rules are considered first.
This is called the Priority, and it has higher precedence than other things which influence rules
  (e.g. biases and weights).
In fact, weights are completely ignored under certain Priorities!

To select a Priority `x`, add it to the beginning of the rules block with `PRIORITIZE(x)`.
Priorities may have extra arguments, e.g. `PRIORITIZE(x, 1, "hello")`,
  but none of the built-in ones actually use this feature.

* `PRIORITIZE(everything)` : **The default priority**.
Each potential rule application has the same chance of being selected,
   at least before weights and biases are applied.
So if rule A has 100 possible applications and rule B has only one,
  then rule A is 100x more likely to be chosen than B.
If rule B has `*2` weighting, then A drops to 50x more likely.
* `PRIORITIZE(fair)` : Every rule has the same chance of being picked,
  at least before weights and biases are applied.
This means that a rule with 200 possible applications and one with 2 possible applications
  would have the same odds of being chosen.
* `PRIORITIZE(earliest)` : Always choose the first rule that has at least one possible application.
Rule weights are ignored under this Priority.
* `PRIORITIZE(latest)` : Opposite of `earliest`, always choose the last rule that has at least one possible application.
Rule weights are ignored under this Priority.
* `PRIORITIZE(common)` : Always choose the rule with the most potential applications.
The counts for each rule are scaled by that rule's weight.
This is similar to the default `everything` behavior, except it *always* chooses the most popular rule.
* `PRIORITIZE(rare)` : Always choose the rule with the fewest potential applications.
The counts for each rule are scaled by that rule's weight.

To define your own priorities, define a new `AbstractMarkovRewitePriority`
  and implement the interface described in its doc-string.
Custom priorities can take extra arguments with the syntax `PRIORITIZE(my_priority, ...)`.
They can also force the rewrite op to end immediately,
  by rule an invalid rule index or a rule that doesn't match anywhere.

### Threshold

By default `@rewrite` will apply its rules forever -- until there are no matches left.
However you can provide a *threshold* as the first argument, to limit this.

* If you want to restrict it to a hard-coded number of matches, pass that number:
  `@rewrite 10 R=>G`.
* If you want to make it relative to the total number of pixels/voxels in the grid,
  pass a simple multiplication or division statement (**in parentheses** if using the simpler macro syntax):
  `@rewrite (area/10) R=>G`. It's automatically rounded and clamped >=1.
* If you want to make it relative to the average length of the grid along each axis,
  pass a simple multiplication or division statement (**in parentheses** if using the simpler macro syntax):
  `@rewrite (0.5*length) R=>G`. It's automatically rounded and clamped >=1.
* If you want a randomized threshold value, pass a range: `@rewrite (area/100):20`.
It's OK if the "max" is greater than the "min".

### Symmetry

By default each rule can be applied along any grid axis, and in either direction
  (`+` meaning the first pixel of the rule is at the min end of the strip; `-` meaning it's at the max end).
To change this, add the following modifier to the end of the rule: `\[ axis1, axis2, -(axis3)... ]`,
  naming each matchable axis (using `x`/`y`/`z`/`w` and `1`/`2`/... interchangeably),
  and optionally a single direction with `+(a)` or `-(a)`.

Due to limitations of Julia's parser, parentheses are **required** when giving a direction for a numeric dimension,
  e.g. `+(2)` vs `+y`.
Otherwise `+2` is parsed as just `2` and permits either direction!

Some examples:

````julia
# Single-rule example
@rewrite RR => Yb  \[-x, Y]  # Only allowed to face -x, -Y, +Y

# Multi-rule example
@rewrite begin
    RR => YY            # May run along any orthogonal direction
    RR => Tb \[+(3)...] # Only able to point upward in the third, fourth, etc axes
    RG => GR \[x, Y]    # Runs along any 2D direction, but e.g. on 3D grids it cannot point along the Z
end
````

Symmetry axes put a lower-bound on the grid's dimensionality.
For example, if you mention a Z axis in any rule's symmetry
  then the whole algorithm cannot run on a 1D or 2D grid.

This modifier must be written **after** the [weight modifier](#bias-and-weights) (see below).

Note that symmetries are automatically restricted to prevent redundant matches.
In particular:
  * If the rewrite rule is a single pixel it will always be given the symmetry `\[ +(x) ]` --
  any symmetry statement you provide is ignored.
  * If the rewrite rule is symmetric along the strip in both source and destination
  (meaning symmetries `-a` and `+a` are equivalent), then all symmetries are converted to the positive direction --
  `\[ -x, y ]` becomes `\[ +x, +y ]`.

### Bias and Weights

As mentioned before, a `@rewrite` with multiple rules will normally choose what to do with uniform randomness.
You can change this behavior by adding weights to each rule,
  and/or a "bias" section to the end of the statement (as either a single expression or a block):

````julia
# A single bias term:
@rewrite R=>G   field(R->Y & B)

# Weighted rules, no bias term:
@rewrite begin
    R => G  *2   # Weight of 2, twice as common
    Y => B       # Default weight of 1
    B => Y  /2   # Half as common; could also write *0.5
end

# Several bias terms:
@rewrite R=>G  begin
    field(R->Y & B)
    field(G->R, live)
end
````

The various biases are described below.

#### Weighted rules

When a `@rewrite` has multiple rules to choose from, then by default each has a weight of 1.
You can change this weight to any positive value.
For example if a weight was changed to `0.2` then it's one-fifth as likely to be chosen
  than under a uniform weighting.

Weights are stated as a multiplication or division right after the rule,
  and only allowed in blocks (since weight has no effect on a single rule anyway):

````julia
@rewrite begin
    # Usually replace G with Y, sometimes replace it with R
    G => Y  *3
    G => R
end
````

Weight modifiers must be specified after the mask (e.g. `%0.5  *3`)
  and before the [symmetry modifier](#symmetry) (e.g. `*3  \[ -x ]`).

#### Bias: Field

> *This is our answer to the original MarkovJunior's `<field>` node.*

A bias of `field(A)` will first compute a **pathing field** outward from all `A` cells,
  storing integer distance from the nearest `A`.
Then all rewrites are biased to stay as low as possible in this field, close to the original `A` cells.
There are many ways to extend this behavior.

* Provide multiple colors for the source cells: `field(RGBwbg)`.
* Flip the bias to be *away* from source cells, using a minus sign: `field(-A)`.
* Force the field to be recomputed any time paths are affected (e.g. when source cells are rewritten)
  with a `live` parameter: `field(A, live)`.
  * Without this parameter the field is computed once when the rewrite op starts,
   which is much more performant.
  Don't recompute unless you really need it!
* Limit the cells in the pathing field by providing a "path": `field(A->BC)` will focus on paths
  from `A` cells through `B` and/or `C` cells;
  others are "outside" the pathing field and cannot be rewritten.
  * If you want outside cells to not be forbidden but merely have the lowest Bias,
  pass a `soft` parameter: `field(A->BC, soft)`.
    * If using randomness (see below), you may want to change how severe the penalty is for outside cells: `field(A->BC, soft=3.5)`. The default is 2.
  * When flipping the bias of a field with paths, like `field(-A->BC)`,
  a friendlier syntax is available: `field(A<-BC)`.
  Note that if a field is reversed, has paths, and is `soft`, then outside cells are rewritten *first*! And the randomization penalty (e.g. `soft=3.5`) becomes a bonus.
  * To flip the pathing so that only writes *outside* the path are allowed, pass `forbidden`:
  `field(A->BC, forbidden)`. Note that in this mode it makes no difference whether the path itself is flipped.
* A field that has paths, like `field(A->BC)`, can also be assigned an "anchor": `field(A->BC & DE)`.
  Cells are outside the pathing field if they don't have a connection to an anchor --
  in this case, to a `D` or `E` cell through `A`, `B`, and `C` cells.
  * With anchors, the distinction between source and path cells becomes important:
  an anchor connection doesn't pass through source cells unless they're also listed as path cells!
* To switch from deterministic Bias to weighted-random, pass `randomness=x` (or `random` which is equivalent to `randomness=0.5`).
  The value ranges from 0 to 1:
  * 0 is "not random", the default behavior.
  * 0.5 is linear weighted-randomness, where an option 3x farther along the path is 3x less likely to be chosen (or *more* likely if the path is flipped).
  * 1 is uniform-randomness, where distance along the path does not matter (and neither does being off the path if `soft` is true).
A value of exactly 1 is equivalent to disabling the bias.
  * When `soft` is true, options off the path are considered to have half as much weight
  as the worst options on the path.
* To help mix this bias with others, pass `scale=x` (default is 1).
Higher values make this bias more important.
* The final bias value for one potential rewrite is an average
of the individual biases for each affected pixel.
You can change the combination logic by passing `combo=X`.
These different modes are *roughly* normalized,
  so the bias always has about the same magnitude range of 0 to max-path-length:
  * `average` (**default**)
  * `max`
  * `min`
  * `first` (the min pixel of the rewrite rule)
  * `deviation` (the amount of variation between pixels in the rewrite)
  * `diff` (max bias - min bias)

#### Bias: Solve

A bias of `solve(b->w)` will try to turn each `b` pixel into `w`.
This is our version of original Markov Junior's "inference" node (a.k.a. `<observe>`).

**TODO: Finish**

#### Bias: Temperature

This bias adds randomness to other biases, where `f` is the scale of the randomness.
More precisely, a uniform-random value from 0 to `f` is added to the total Bias for each rewrite choice.

Temperature works great for small-scale randomness, shuffling the most preferred moves
  with the almost-most-preferred moves.
It does not help for large-scale randomness, because low-bias moves cannot be chosen
  except by raising the temperature so high that the other biases are meaningless.
Instead, individual biases need to offer their own way of doing large-scale randomness.

### Multidimensional matches

A tweak to the syntax for rewrite rules allows you to match *blocks* of pixels -- 2D, 3D, and more!
Instead of a flat string of characters like `a[bc]d_e[fgh]`, use Julia's multidimensional arrays, e.g.

````julia
# A 2D block rule:
[ a b c
  d e f
] => [
  b a f
  e d c
]

# A 3D block rule, using the usual rewrite tricks:
[ a   [ab]    c
  d    e      f ;;;
  g    h      i
  j    k      l
] => [
  c    _      a
  d    e      f ;;;
  g [2, 1, 1] i  # Now taking by index requires 3 indices -- X, Y, Z!
  j  {jkl}    l
]

# A 4D block rule, with modifiers:
[ a b c
  d e f ;;;
  g h i
  j k l ;;;;
  m n o
  p q r ;;;
  s t u
  v w x
] => [
  c b a
  e d f ;;;
  g g i
  l k l ;;;;
  m n o
  p q r ;;;
  s t u
  v w x
]  *2  \[
    # Symmetry is now per-block-axis; read the section about it further below.

    #   Horizontal (X and Y) axes of the block must stay horizontal,
    #     but can rotate and flip amongst each other.
    x[ x, y ],
    y[ x, y ],

    #   Z axis of the block must stay pointed upwards.
    z[ +Z ]

    #   Fourth axis of the block must stay along the fourth axis of the grid,
    #      but may flip backwards along it.
   # w[ W ]
    # Does not actually need to be written because it's implicitly the only choice left for a 4D grid.
]
````

> *If you're familiar with Julia, be aware that we flip the first two array axes internally.*
> *Mathematicians think in terms of Row->Column, but we think about X->Y, so for us Column should be the first axis.*

Blocks can have fewer dimensions than the grid itself, but not more.
For example a 3D block can be used on a 4D grid but not a 2D grid.

### Multidimensional symmetry

By default, like with normal rewrite rules, blocks can be oriented any which way.
This gets very complicated as you go beyond 2D, but there is a simple way to think about it:
  a block's orientation on the grid is merely a process of selection.

1. Pick a grid axis and sign (towards + or towards -), for the first block axis to orient along.
2. Pick a *different* grid axis and a sign, for the second block axis to orient along.
3. Repeat for the other block axes.

This simple model leads to a surprisingly simple description of multidimensional symmetries:
For each block axis, provide an explicit list of the axes/directions it may choose to orient along.
The axis choices are made in the order you write,
  and afterwards any unmentioned block axes can choose any available grid orientation.

Here are examples of block symmetry modifiers:

````julia
# For a 2D block that only allows flipping of the X axis:
\[   x[ x ], y[ +y ]   ]
# If we assume the grid is always 2D, it can be simplified:
\[   y[ +y ]   ]

# For a 3D block that allows horizontal rotation/flips but no messing with the vertical:
\[ x[x, y], y[x, y], z[ +z ] ]
# Again, If we assume the grid is always 3D then it shortens a lot:
\[z[+z]]
````

> **NOTE**: *you must add commas between the elements of a symmetry statement!*
> *Otherwise the parser will break and throw a strange parsing error,*
> * which is hard for us to catch and rethrow more clearly.*

As with the simpler "strip" rewrite rules, you can use `a...`
  to indicate matching with any extra dimensions starting at `a`.

````julia
\[
    # Block's Y axis can face anywhere but the grid's Y axis.
    y[ x, z... ]
]
````

#### Grouping axes

When two or more axes have a close relationship, you may group them together
  and provide explicit permutations as tuples.
For example, this permits the block X and Y axes to rotate amongst themselves but never flip.

````julia
\[  (x, y)[ (+x, +y), (+y, -x), (-x, -y), (-y, +x) ]  ]
````

> *This use-case actually has its own macro, `\[  (x, y)[ @rotations(x, y) ]  ]`.*
> *The macro can be used with any number of axes too!*

This next example either allows the block's X to flip,
  or for the whole block to reorient along the grid's WZ plane.

````julia
\[
    (x, y)[ (x, +y), (z, w), (w, z) ]
]
````

Here is a 3D block symmetry that only allows full inversion of the block:

````julia
\[
    (x, y, z)[ (+x, +y, +z), (-x, -y, -z) ]
]
````

#### No-flipping axes

One last feature is the ability to forbid sets of block axes from containing a flip,
  no matter where they end up, with the syntax `{axes...}`.
For example if a 2D block should be able to orient anywhere in a grid but only through rotations,
  you can write `\[ {x,y} ]`.

Suppose we're generating a 4D grid where the first axis is meant to be Time.
The below 3D rewrite rule can be rotated any way in space --
  not flipped, and none of its axes are allowed to stretch through Time.
````julia
[ a b c
  d e f ;;;
  g h i
  j k l
] => [
  c b a
  d e f ;;;
  g h i
  j k l
] \[
    {x, y, z},
    x[ y, z, w ],
    y[ 2, 3, 4 ],
    z[ 2, z, 4 ],
]
````

#### Symmetry macros

To simplify the most common use-cases, some macros are provided that inject permutations.

* `@rotations(axes...)` inserts the rotation permutations of two or more axes.
You can limit the first two axes to only rotate between themselves by doing `\[   (x, y)[ @rotations(x, y) ]   ]`.
You can limit the block axes `y`, `z`, and `w` to only rotate between three grid axes `x`, `y`, and `z`
  by doing `\[   (y, z, w)[ @rotations(x, y, z) ]     ]`.
* `@inversions(axes...)` inserts permutations where one or more axes flip.
For example `\[  (x, y)[ (+y, +z), @inversions(y, z) ]   ]` turns into
`\[   (x, y)[ (+y, +z), (-y, +z), (+y, -z), (-y, -z) ]   ]`.

## `@fill`

To make a deterministic modification to a grid, you can use
  `@fill 'C' [S(A=N, B=M)] [rule] [mask]`.

This operation will fill a single color within a specific region of the grid.
If the region (`S(...)` above) isn't provided, it defaults to the entire grid.
A rule can optionally filter which color pixels are affected,
   so you can also think of it as a "replace" operation.

If this box's space is described with vectors, and not scalars,
  then the box can't be higher-dimensional than the grid!
However it can be lower-dimensional -- for additional grid dimensions the box stretches along the entire space.

The parameters are as follows:
* `'C'` is the color to draw.
* `S` is the kind of space covered by the box, either:
  * `uv` for continuous space where 0 is the grid's min corner, 1 is the max corner,
and box size is always at least one pixel (so you can pass `size=0` to draw exactly one pixel).
  * `pixel` for integer pixel coordinates (starting at 1).
You can still provide fractional values, in which case coordinates will be rounded.
* `A` and `B` are the parameters defining the box, chosen among the following:
  * `min`
  * `max` (inclusive)
  * `size`
  * `center`
* `X` and `Y` are the values for those parameters, either:
  * A single number, representing the value along all axes
  * A vector like `(3, 20)`, *which also forces the entire grid to have at least that many dimensions*.
* `[rule]` is an optional condition for which pixels are affected, using `+` to whitelist and `-` to blacklist.
For example `-RGB` means to affect every pixel in the box except Red Green and Blue;
`+RGB` means to affect *no* pixels in the box except Red Green and Blue.
* `[mask]` is an optional mask statement which randomly forbids some of the pixels in the same way as for [rewrite rules](#), for example `%0.75` forbids 25% of all pixels.

By default this op will sometimes finish in a single tick, and sometimes spread itself out over several ticks.
You can force it to always complete in a single tick by adding a pragma to the algorithm: `@pragma fast_fills`.

Note that a 1D vector like `(5, )` has different behavior than a scalar like `5`
  when extrapolating to higher dimensions.
The 1D vector lets every axis from 2 onward take up the entire grid,
  while the scalar value is used for every axis.

## `@sequence`

````julia
@sequence (length/10) begin
    # Redden 20% of the pixels in the top-left corner of the grid.
    @fill 'R' uv(min=0, size=0.5) %0.2
    # Turn one of them into green.
    @rewrite 1 R => G
    # Mark a direction pointing outward from the green pixel.
    @rewrite G_ => GB
end path(G<-R<-w) # Add a pathing bias
````

The op `@sequence [threshold] begin ... end` describes a chronological sequence of actions,
  similar to the algorithm itself but with Thresholding options.
For example:

The `threshold` is identical to [the Threshold for `@rewrite` statements](#threshold),
  determining how many times to run, but with an extra option:
 `repeat` makes the sequence repeat until the first inner operation fails to do anything.

The last argument to a sequence is a bias statement or block of statements,
  identical to [bias for `@rewrite` statements](#bias-and-weights).
These biases are inherited by all operations within it
  (primarily `@rewrite` but other ops may respect `bias` settings too).

## `@upscale`

````julia
# Simply convert each pixel into a 2x2 block of itself
@upscale (2, 2)

# Convert each pixel into a 3x2 block according to the following rules
@upscale (3, 2) begin
    R => [
        w b g
        M M M
    ]
    G => [ M;M;M ;; Y;Y;Y ] # Alternative Julia syntax for multidimensional arrays
    B => [
        # Each element that names multiple colors will pick a random one each time
        wb wb wb
        g  g  g
    ]
    # All other pixels become a 3x2 block of themselves
end
````

This increases the size of the grid by stretching each pixel a certain number of times along each axis.
You may also add rewrite rules that replace specific pixel colors with specific patterns.

The basic syntax is `@upscale (factors...) [rewrite rules]`.
Factors are the stretch amounts along each axis, for example `(2, 2)`.

The stretch factor for any extra grid axes defaults to 1, a.k.a. no change.
However if you add an ellipsis to the end, then those extra axes all take on the last factor.
For example `@upscale (2...)` will double the size along *all* axes.

The rewrite rules can be nothing, or a single rule, or a `begin ... end` block of rules.
If no rules match a source pixel, that pixel keeps its value during upscaling.
Rules are formatted as a source color, arrow, and destination block,
  with individual destination values being a single color OR a group whose value is picked randomly each time.
For example:

````julia
# This rule only works in a grid of at least two dimensions,
#   and a rewrite rule which uses a factor of 2 for the first two axes.
# For every grid dimension above 2D, the block is extruded along that dimension.
# We borrow the 'set' syntax from rewrite rules to allow random outcomes for certain pixels.
G => [
  R  BG
  G   M
]
````

If the grid upscales along more dimensions than the block, that block is extruded along the extra dimensions.
For example, in the above example rule, a 3D upscale would use that same 2D slice at each Z level.

## `@downscale`

````julia
# Simply replace each 2x2 block of pixels with a single pixel,
#   chosen randomly from that block.
@downscale (2, 2)

# Replace each 3x2 block of pixels with a single pixel,
#   usually with random choice but certain blocks have more specific outcomes.
@downscale (3, 2) begin
  # A row of RGB, followed by a greyscale row bookended with white and black,
  #   becomes a Red pixel.
  R <= [
    R   G   B
    w  wgb b
  ]
  # A block of Magenta becomes either Green or Blue.
  GB <= [ M;M;M ;; M;M;M ]

  # Any other blocks fall back to the usual behavior of picking a random pixel.
  # If a grid has a size not evenly divisible by our factor (3x2),
  #    there are partial blocks at its ends which all get the default (random) behavior.
end
````

This decreases the size of the grid by splitting it into regular blocks of pixels,
  and shrinking each block to a single pixel.
By default the shrunken pixel chooses a random color from its block,
  but you may add rewrite rules to catch certain patterns and map to a specific color.

The basic syntax is `@downscale (factors...) [rewrite rules]`.
Factors are the box sizes (i.e. down-scale amounts) along each axis.

The factor for any unmentioned/extra grid axes defaults to 1, a.k.a. no change.
However if you add an ellipsis to the end of your list, then those extra axes all take on the last factor.
For example `@downscale (2...)` will halve the size along *all* axes.

The rule syntax is destination, then left-arrow, then source block.
The source block must match the factor (e.g. downscaling by `(3,2)` requires that each rule have a 3x2 block),
  though extra dimensions can be omitted and your block will simply be extruded along that axis.

Both destination pixel and source pixels can use Set syntax,
  listing multiple colors at once to indicate a random choice among them.

## `@convolve`

Convolution is the technical term for processing an image
  by replacing each pixel with a particular blend of itself with its neighbors.
The specific blend being used is called the Kernel.
For example, blurs are defined this way and you get different kinds of blurs by using different Kernels.
It's also a popular way to implement Cellular Automata such as the Game of Life.

With MarkovJunior's extremely limited palette, convolution is defined in a simpler way:

* The Kernel is an array of which neighboring pixels to include, optionally only including particular colors at particular neighbor points.
* Different kernels can be used depending on the color of the central pixel. One color can even define multiple kernels, chosen randomly per pixel.
* For each pixel, a Kernel is chosen and all included neighbors are counted by their color. 
* A series of if-then statements is used to pick the new value based on integer math with the color counts (and with the color of the central pixel).

**TODO: syntax**

## `@not_animated`

This is a small helper to tweak rendered animations.

During normal algorithm runs, for optimization purposes, certain complex ops (`@rewrite`)
  may change the order in which things happen without affecting the final result.
However when rendering an algorithm to video,
  these optimizations are disabled so that the frame-by-frame behavior looks correct.

The nested operator `@not_animated [op]` re-enables those optmizations
  for the duration of the given op.