# MarkovJunior.jl

A Julia reimagining of [this awesome procedural generation algorithm](https://github.com/mxgmn/MarkovJunior/),
  able to generate in any number of dimensions.
It is planned for release as both a standalone executable to create scenes/images/video,
  and a C-like library that can be used anywhere.

````julia
# In the Julia REPL:
] add MarkovJunior
using MarkovJunior
markovjunior_run_tool()
````

![A screenshot of the Terrain.jl scene](docs/Screenshot.png)

[![A demo video, sped up 2x](https://img.youtube.com/vi/_I9aU8Ad-tE/0.jpg)](https://www.youtube.com/watch?v=_I9aU8Ad-tE)

*A YouTube video demo*

## Description

The released 0.1 version was a proof-of-concept, with a 2D renderer --
  totally functional and fun to play with.

The current version has lots of nice features:

* A terse, extendible, broadly-featured DSL
* Both 2D and 3D rendering
* More kinds of Operations and Biases

[You can find documentation on the new syntax here](docs/dsl.md).

Math and rendering is all built on top of my [B+ game framework](https://github.com/heyx3/B-plus).

The unit tests in *test/runtests.jl* are mostly focused on parsing,
  while correctness is mostly checked by running sample scenes in the GUI,
  but there are exceptions for important algorithm-heavy internal functions.

## Why higher-dimensional grids?

I decided to support any-dimensional grids, mainly because it was cool and Julia makes it feasible.
However there are real, if underexplored, benefits!

* You can add a Time axis to the grid to get animated output
* You can add an extra axis to store multiple pixels of data at each output pixel, like temp variables.
* You could perhaps use it to foster logical connections between disparate areas of the grid,
  by connecting them through the extra dimensions.

If you come up with any higher-dimensional tricks, share them in a Github Issue!

## Development

On the main branch, this tool uses the main branch of B+ (and its sub-packages).
Clone BplusCore, BplusApp, BplusTools, and Bplus,
  then run the following from the Julia REPL:

````julia
# Add local B+ sub-packages to B+
] activate ../Bplus.jl
] dev ../BplusCore.jl ../BplusApp.jl ../BplusTools.jl
# Add local B+ to MarkovJunior
] activate .
] dev ../Bplus.jl
````

Due to a [bug with FixedPointNumbers.jl](https://github.com/JuliaMath/FixedPointNumbers.jl/issues/317)
  -- possibly a [bug with PackageCompiler itself](https://github.com/JuliaLang/PackageCompiler.jl/issues/1092) --
  the public release can't be used with PackageCompiler.
Instead, you need to [clone this fork](github.com/Kyjor/FixedPointNumbers.jl/tree/bugfix/remove-asserts-from-precompile),
  at branch `bugfix/remove-asserts-from-precompile` (e.g. `git checkout bugfix/remove-asserts-from-precompile`),
  and link it to this project:

````julia
# Add it to the local Julia environment:
] activate
] dev ../FixedPointNumbers.jl
# Add it to this project:
] activate .
] dev ../FixedPointNumbers.jj
# Not sure why both of the above are needed, but they are.
````

## Scenes

Different algorithm setups can be found in the *scenes/* folder.
They are commented with explanations of how they work.

A small handful of them are still not converted to the v0.2 syntax,
  so don't be alarmed if they fail to run!