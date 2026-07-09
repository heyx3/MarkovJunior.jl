# v0.3

* Add 3D view to the scene renderer
* Support compilation to a GUI-tool exe, DLL, and IPC service exe.
See *scripts/compile_standalone.jl*
* Add more features to existing Ops
* Add `field()` Bias (see [DSL documentation](docs/dsl.md#bias-field))
* Fix up the interface for custom extensions (see [*add_ons* documentation](docs/add_ons.md))
* Add "reused" allocator and use it as the default
* GUI improvements

# v0.2

* Switch to new DSL
* Redesign the code architecture to go with the new DSL and better fit the underlying algorithm
* Update the built-in scene files

# v0.1

Initial version, proof-of-concept syntax and algorithm features.