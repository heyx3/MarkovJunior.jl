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

Lots of nice features have already been implemented:

* A terse yet broadly-featured DSL
* Both 2D and 3D rendering of scenes
* Any-dimensional rewrite rules with complex symmetry specifications
* Numerous hooks to extend the DSL -- new Ops, Biases, and more
* Building to a standalone exe and use-anywhere DLL library

[You can find documentation on the new syntax here](docs/dsl.md).

Math, GUI, and rendering are all built on top of my [B+ game framework](https://github.com/heyx3/B-plus).

## Why higher-dimensional grids?

I decided to support any-dimensional grids, mainly because it was cool and Julia makes it feasible.
However there are real, if underexplored, benefits!

* You can add a Time axis to the grid to get animated output.
* You can add an extra axis to store multiple pixels of data behind each output pixel, like temp variables.
* You could perhaps use it to foster logical connections between disparate areas of the grid,
  by connecting them through the extra dimensions.

If you have any thoughts on higher-dimensional tricks, share them in a Github Issue!

## Scenes

Different algorithm setups can be found in the *scenes/* folder.
They are heavily commented to help you learn.

A small handful of them are still not converted to the v0.2 syntax,
  so don't be alarmed if they fail to run!

## Build size

Due to Julia's unique JIT architecture, it essentially needs a whole compiler inside its runtime.
This means that the executabls and dll's tend to be much larger than is ideal,
  and also do not play well with mobile.

Fortunately it is possible to get around these problems once you decide you need to,
  by building your own Julia package which has ours as a dependency.

* Cut out the GUI tool and its dependencies to greatly simplify the dependency tree.
* Pre-parse all the MarkovAlgorithm instances you actually need to use, store them in your package,
  and make sure they all get precompiled.
* Make the build in a similar way to how this package does it, but with some parameter changes:
  * Pass `filter_stdlibs=true`. If there are any runtime crashes due to missing libs,
  add those libs as a direct dependency to your project.

This package has plans to reduce binary size in the future.

## Development

On the main branch, this package uses the main branch of B+ (and its sub-packages).
So you should clone BplusCore, BplusApp, BplusTools, and Bplus,
  then run the following from the Julia REPL:

````julia
# Add local B+ sub-packages to B+
] activate ../Bplus.jl
] dev ../BplusCore.jl ../BplusApp.jl ../BplusTools.jl
# Add local B+ to MarkovJunior
] activate .
] dev ../Bplus.jl
````

### Compilation (to .exe or .dll)

We use the standard Jula project *PackageCompiler.jl* to compile the main codebase into a DLL,
  and the codebase+GUI tool into an EXE.
For convenience, *scripts/compile_standalone.jl* encapsulates all the project-specific logic.
You must either pass `-exe` or `-dll` to make it build an executable or a library, respectively.

### FixedPointNumbers.jl workaround

Due to a [bug with FixedPointNumbers.jl](https://github.com/JuliaMath/FixedPointNumbers.jl/issues/317)
  -- possibly a [bug with PackageCompiler itself](https://github.com/JuliaLang/PackageCompiler.jl/issues/1092) --
  the public version of that package can't be used with PackageCompiler.
Instead, you need to [clone this fork](github.com/Kyjor/FixedPointNumbers.jl/tree/bugfix/remove-asserts-from-precompile),
  at branch `bugfix/remove-asserts-from-precompile` (e.g. `git checkout bugfix/remove-asserts-from-precompile`),
  and link it to this project:

````julia
# First add it to your global Julia environment:
] activate
] dev ../FixedPointNumbers.jl
# Second add it to this project:
] activate .
] dev ../FixedPointNumbers.jj
# Not sure why both of the above are needed, but they are.
````

### `nvpatch` workaround

In the GUI tool, due to the use of *Bplus.jl*, we need to run on discrete GPU's instead of integrated ones.
The best way to ensure graphics drivers know this is to export specific constants in the compiled executable.
Unfortunately there's no way to do that through *PackageCompiler.jl*, so we use an external tool
  [called `nvpatch`](https://github.com/toptensoftware/nvpatch).

If you build this project into an executable without `nvpatch` installed, you'll get a stern warning.
Users of the GUI tool will have to force the discrete GPU themselves
    through Nvidia's Control Panel (or AMD's equivalent) to avoid errors on tool startup.

### Testing

The biggest test in *test/runtests.jl* is focused on parsing,
  while correctness is mostly checked by running all sample scenes in the GUI.