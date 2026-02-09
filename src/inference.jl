########################
#  Path inference

"
A grid constraint based on some cells' ability
  to reach other cells, through some third set of cells.

Biases the 'source' types to use rule applications
   that push through the 'path' types toward the 'dest' types.
If `invert` is true then they're biased to push away instead.
"
struct InferPath
    source_types::CellTypeSet
    dest_types::CellTypeSet
    path_types::CellTypeSet
    invert::Bool
    temperature::Float32
    recompute_each_time::Bool

    InferPath(src, dest, path;
              invert::Bool = false,
              temperature::Float32 = 0.0f0,
              recompute_each_time::Bool = true) = new(
        if src isa CellTypeSet
            src
        else
            CellTypeSet(src)
        end,
        if dest isa CellTypeSet
            dest
        else
            CellTypeSet(dest)
        end,
        if path isa CellTypeSet
            path
        else
            CellTypeSet(path)
        end,
        invert, temperature, recompute_each_time
    )
end

dsl_string(p::InferPath) = string(
    "@path ", @sprintf("%f", p.temperature), " ",
    p.recompute_each_time ? "recompute " : "",
    p.invert ? "invert " : "",
    '"', dsl_string(p.source_types), "\" => \"",
         dsl_string(p.path_types), "\" => \"",
         dsl_string(p.dest_types), "\" "
)

"The current state of an `InferPath` constraint applied to a specific grid"
struct InferPath_State{N}
    constraint::InferPath
    grid::CellGrid{N}
    potential::Array{UInt32, N} # Max possible value represents 'unreachable'

    search_frontier_buffer::Dict{CellIdx{N}, UInt32}

    function InferPath_State(constraint::InferPath, grid::CellGrid{N}) where {N}
        potential = fill(typemax(UInt32), size(grid))

        output = new{N}(constraint, grid, potential, Dict{CellIdx{N}, UInt32}())
        recalculate!(output)

        return output
    end
end

"Regenerates this InferPath constraint's `potential` field"
function recalculate!(ips::InferPath_State{N}) where {N}
    @markovjunior_assert(isempty(ips.search_frontier_buffer))

    LOG_INFERENCE_PATH && println(
        "Calculating ", dsl_string(ips.constraint), "...\n",
        "\tSource=", collect(ips.constraint.source_types),
          " | Paths=", collect(ips.constraint.path_types),
          " | Dest=", collect(ips.constraint.dest_types))

    # Seed the potential field by setting the dest cells' potentials to 0.
    fill!(ips.potential, typemax(UInt32))
    for grid_idx::CellIdx{N} in one(CellIdx{N}):vsize(ips.grid)
        if ips.grid[grid_idx] in ips.constraint.dest_types
            ips.search_frontier_buffer[grid_idx] = UInt32(0)
            LOG_INFERENCE_PATH && println("\tDest at ", grid_idx)
        end
    end

    # Apply the search frontier to the potential field,
    #    each time queueing up neighboring cells to be re-evaluated.
    while !isempty(ips.search_frontier_buffer)
        (grid_idx, new_value) = pop!(ips.search_frontier_buffer)
        LOG_INFERENCE_PATH && println(
            "\tProcessing ", grid_idx,
            " | previous=", ips.potential[grid_idx], " | new=", new_value
        )

        @markovjunior_assert(new_value < ips.potential[grid_idx])
        ips.potential[grid_idx] = new_value

        neighbor_new_value = new_value + UInt32(1)
        LOG_INFERENCE_PATH && println(
            "\t\tNeighbors could have potential ", neighbor_new_value
        )
        for neighbor_axis in 1:N
            for neighbor_dir in -1:2:1
                neighbor_axis_pos = grid_idx[neighbor_axis] + neighbor_dir

                if (neighbor_axis_pos > 0) && (neighbor_axis_pos <= size(ips.grid, neighbor_axis))
                    neighbor_grid_idx = let i = grid_idx
                        @set i[neighbor_axis] = neighbor_axis_pos
                    end
                    LOG_INFERENCE_PATH && println(
                        "\t\tChecking neighbor at ", neighbor_grid_idx, " | value=", ips.grid[neighbor_grid_idx],
                            " | officialPotential=", ips.potential[neighbor_grid_idx],
                            " | upcomingPotential=", get(ips.search_frontier_buffer, neighbor_grid_idx, NaN),
                            " | ", ips.grid[neighbor_grid_idx] in ips.constraint.path_types,
                                   "/", neighbor_new_value < get(ips.search_frontier_buffer, neighbor_grid_idx, ips.potential[neighbor_grid_idx])
                    )
                    if (ips.grid[neighbor_grid_idx] in ips.constraint.path_types) &&
                        (neighbor_new_value < get(ips.search_frontier_buffer, neighbor_grid_idx, ips.potential[neighbor_grid_idx]))
                    #begin
                        ips.search_frontier_buffer[neighbor_grid_idx] = neighbor_new_value
                        LOG_INFERENCE_PATH && println(
                            "\t\t\tOur new potential is better! Adding it to the search frontier!"
                        )
                    end
                end
            end
        end
    end

    LOG_INFERENCE_PATH && begin
        println("Final potential:")
        display(ips.potential') # Transposed to match the display
        println("\n")
    end

    return nothing
end

"Gets the delta-potential of a given rule application according to this path"
function infer_weight(inference::InferPath_State{N},
                      rule::CellRule, at::CellLine{N},
                      cached_current_types_along_line::Vector{UInt8},
                      cached_current_types_in_line::CellTypeSet,
                      rng::PRNG
                     )::Float32 where {N}
    delta_potential::Float32 = rand(rng, Float32) * inference.constraint.temperature
    for_each_cell_in_line(at) do idx::Int32, pos::CellIdx{N}
        was_path_src::Bool = inference.grid[pos] in inference.constraint.source_types
        is_path_src::Bool = exists(rule.output[idx]) && (rule.output[idx] in inference.constraint.source_types)
        if was_path_src || is_path_src
            raw_potential = inference.potential[pos]
            potential = if raw_potential == typemax(UInt32)
                -1.0f0
            else
                convert(Float32, raw_potential)
            end
            delta_potential += potential *
                                 (is_path_src - was_path_src) *
                                 (inference.constraint.invert ? -1.0f0 : 1.0f0)
        end
    end
    return delta_potential
end


########################
#  Overall inference

"
Various high-level constraints on a grid during generation:
 * *Paths* bias certain colors ('source') to be placed towards certain other colors ('dest')
along paths made by a third set of colors ('path').
 * *Temperature* reduces the strength of inference by randomizing its outputs.

When constructing an instance by concatenating several other ones,
  the last instance is used for unique fields like `temperature`.
"
struct AllInference
    paths::Vector{InferPath}
    paths_by_type::Dict{UInt8, Vector{Int}} # Each value element is an index into 'paths'

    temperature::Float32

    function AllInference(paths::Vector{InferPath}, temperature::Float32)
        paths_by_type = Dict{UInt8, Vector{Int}}()
        for (i, path) in enumerate(paths)
            for type in path.source_types
                push!(
                    get!(() -> Vector{Int}(), paths_by_type, type),
                    i
                )
            end
        end

        return new(paths, paths_by_type, temperature)
    end
    function AllInference(to_merge::AllInference...)
        paths = InferPath[ ]
        paths_by_type = Dict{UInt8, Vector{Int}}()
        temp = 0.0f0

        path_idx_offset = 0
        for tm in to_merge
            append!(paths, tm.paths)
            temp = tm.temperature

            for (type, path_idcs) in tm.paths_by_type
                append!(
                    get!(() -> Vector{Int}(), paths_by_type, type),
                    path_idcs .+ path_idx_offset
                )
            end

            path_idx_offset += length(tm.paths)
        end
        return new(paths, paths_by_type, temp)
    end
end

dsl_string(i::AllInference) = string(
    "@infer begin\n",
    "\t", @sprintf("%f", i.temperature), "\n",
    "\t", iter_join(dsl_string.(i.paths), "\n\t")...,
    "\nend"
)

"Checks whether any inference should be happening"
function inference_exists(constraints::AllInference)::Bool
    return !isempty(constraints.paths)
end


"The current state of all constraints for a specific grid"
struct AllInference_State{N}
    source::AllInference
    grid::CellGrid{N}
    paths::Vector{InferPath_State{N}}
    buffer1::Vector{UInt8}

    function AllInference_State(constraints::AllInference, grid::CellGrid{N}) where {N}
        return new{N}(
            constraints, grid,
            InferPath_State.(constraints.paths, Ref(grid)),
            Vector{UInt8}()
        )
    end
end

"
Get the effect on inference of applying the given rule to the given location.

It's recommended to check `inference_exists()`
   and have a fast path in case you don't need to calculate inference at all.
"
function infer_weight(inference::AllInference_State{N},
                      rule::CellRule, at::CellLine{N},
                      rng::PRNG
                     )::Float32 where {N}
    # Cache the colors along the cell line, for convenience.
    empty!(inference.buffer1)
    types_in_line = let r = Ref(CellTypeSet())
        for_each_cell_in_line(at) do idx, cell::CellIdx{N}
            cell_value = inference.grid[cell]
            push!(inference.buffer1, cell_value)
            r[] = push(r[], cell_value)
        end
        r[]
    end

    # Accumulate inference weights from the various constraints.
    weight::Float32 = rand(rng, Float32) * inference.source.temperature
    for path in inference.paths
        weight += infer_weight(path, rule, at,
                               inference.buffer1, types_in_line,
                               rng)
    end

    return weight
end
"
Checks the total potential of a grid cell (mainly for debugging/visualization).
May have a negative value (if the cell's color is off all inference paths).
"
function infer_weight(inference::AllInference_State{N}, at::CellIdx{N}) where {N}
    weight::Float32 = 0
    for path in inference.paths
        weight += let u = path.potential[at]
            if u == typemax(UInt32)
                -1.0f0
            else
                convert(Float32, u)
            end
        end
    end
    return weight
end

"Recalculates any inference constraints which are marked for frequent recalculation"
function recalculate!(inference::AllInference_State)
    for path in inference.paths
        if path.constraint.recompute_each_time
            recalculate!(path)
        end
    end
    return nothing
end