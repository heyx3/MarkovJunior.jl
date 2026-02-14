# MarkovJunior.jl

A Julia reimagining of [this awesome procedural generation algorithm](https://github.com/mxgmn/MarkovJunior/).

![A screenshot of the Terrain.jl scene](docs/Screenshot.png)

[![A demo video, sped up 2x](https://img.youtube.com/vi/mjBu7Omch1s/0.jpg)](https://youtu.be/mjBu7Omch1s)

*A YouTube video demo*

## Description

The version in main is a proof-of-concept, totally functional and fun to play with.
Currently I am working on a much more comprehensive version with:
* New and improved syntax
* Very robust support for any number of dimensions (e.g. 4D for animated 3D structures)
* More features that the original had (though probably never all of them)
* Several cool tricks that the original doesn't have
* The ability to run as both a standalone tool and a C-style DLL.

The math and rendering is all on top of my [B+ game framework](https://github.com/heyx3/B-plus).
This version of the tool is currently written for 2D scenes, with 3D on the horizon,
  but the algorithm (and DSL) can operate in any number of dimensions!

Note that currently this tool uses the master branch of B+ (and its sub-packages),
  so if you're getting compile errors related to B+ then grab it (and its sub-packages)
  and add to this project with `dev ../Bplus`.

## Syntax

I've developed my own Domain-Specific Language ("DSL"), separate from the original MarkovJunior.

The current version is undocumented but pretty simple to work out from the sample scenes.

The new planned version [is documented here](docs/dsl.md).
