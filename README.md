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

## Usage

A specific algorithm is represented by the [macro `@markovjunior`](docs/dsl.md).
The output of this macro is a `MarkovAlgorithm`.

If you have the macro as a string, call `markov_algo_parse(str)::MarkovAlgorithm`.
Convert the algorithm back to a string with `markov_algo_to_string(algo)`,
   but be aware that cosmetic details get lost in the translation.

To start running an algorithm on a new grid, call `state = markov_algo_start(algo, (32, 32), 0xababcdef)`.
The second argument is the grid size (as many dimensions as you want), and the third is the RNG seed(s).

Iterate on the running algorithm with `markov_algo_step(algo, state)`.
You can pass an iteration count as the third parameter to batch-run multiple steps.
You can also call `markov_algo_finish(algo, state)` to immediately run to completion.

Check on a running algorithm's state with `markov_algo_is_finished(algo, state)`.
Get the grid being operated on with `markov_algo_grid(state)`.

## Scenes

Different algorithm setups can be found in the *scenes/* folder.
They are heavily commented to help you learn.

A small handful of them are still not converted to the v0.2 syntax,
  so don't be alarmed if they fail to run!

## Optimizing the standalone builds

**This section is for anyone who has problems integrating our released executable or DLL into a project.**

Due to Julia's unique JIT architecture, it essentially needs a whole compiler inside its runtime.
Due to our GUI tool, there are many sub-sub-dependencies compiled into the package.
As a result the executable and dll are much larger than is ideal and also don't play well with mobile.
Addtionally if you have any other Julia packages that you'd like to use,
  you'd need to manage multiple copies of julia across multiple libraries!

The way to get around these problems is to make a new Julia package,
  taking ours and any others you want as dependencies,
  then building a new library in the same way we build ours.

Such a package could also precompile particular use-cases for your project,
  like pre-parsing all the specific algorithm instances you plan to run.
This removes the JIT overhead of running each MarkovJunior algorithm for the first time,
  and *that* opens up the door to a static only-AoT-compiled build
  which is mobile-friendly!
Unfortunately AoT Julia builds are still uncommon and a bit of an open problem AFAIK.

If you don't need the GUI tool, cut it and its dependencies out to greatly simplify the dependency tree.
In particular the `Bplus` dependency could be reduced to just `BplusCore`,
  with one or two changes to our package source to accomodate.
We have future plans to automate this.

When building your library, I recommend passing `filter_stdlibs=true` in to PackageCompiler.jl
  to shrink things even further; just look out for runtime crashes due to missing libs.
They can be fixed by adding those libs as a direct dependency.

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

### Linking with the dll: `GenLibFromDll`

PackageCompiler can generate a DLL, [but not the .lib file](https://github.com/JuliaLang/PackageCompiler.jl/issues/687) needed to link it to our headers!
Fortunately it is possible to *generate* a .lib by examining what's exposed in the dll and working backwards.
So, after compilation we use a tool which does exactly this: [GenLibFromDll](https://github.com/KHeresy/GenLibFromDll).
A specific version is directly embedded in this repo under the *scripts/* folder, to prevent link rot.

`GenLibFromDll` depends on Visual Studio Build Tools, so you need those installed.
If you don't have them (or are building for a different OS) you'll get a descriptive error message,
  warning you that there is no easy way for users to link with the DLL.

### Testing

The biggest test in *test/runtests.jl* is focused on parsing,
  while correctness is mostly checked by running all sample scenes in the GUI.