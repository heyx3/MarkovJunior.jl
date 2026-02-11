########################
#  Path inference

"
A grid constraint based on some cells' ability to reach other cells, through some third set of cells.
To wit, if a particular color has a path inference put on it:
  1. No rule is allowed to change cells of that color unless the path can reach that cell
  2. Rules are more likely to be applied when they're further along the path (closer to `dest_types`).
You can flip this by passing `invert=true`.

In cases where this path is not available but other ones are for the same cell,
   this path can apply a penalty to lower the chances of applying a rule to this cell.
The penalty should normally be slightly more than 1 (or equal to it to apply no penalty).
"
struct InferPath
    source_types::CellTypeSet
    dest_types::CellTypeSet
    path_types::CellTypeSet
    invert::Bool
    temperature::Float32
    recompute_each_time::Bool
    penalty::Float32

    InferPath(src, dest, path;
              invert::Bool = false,
              temperature::Float32 = 0.0f0,
              recompute_each_time::Bool = true,
              penalty::Float32 = 1.0f0) = new(
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
        invert, temperature, recompute_each_time, penalty
    )
end

dsl_string(p::InferPath) = string(
    "@path ", @sprintf("%f", p.temperature), " ",
    p.recompute_each_time ? "recompute " : "",
    isone(p.penalty) ? "" : "penalty=$(p.penalty)",
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

    search_frontier_new_potentials::Array{UInt32, N}
    search_frontier_buffer::Queue{CellIdx{N}}

    function InferPath_State(constraint::InferPath, grid::CellGrid{N}) where {N}
        potential = fill(typemax(UInt32), size(grid))

        output = new{N}(constraint, grid, potential,
                        copy(potential), Queue{CellIdx{N}}())
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
    fill!(ips.search_frontier_new_potentials, typemax(UInt32))
    for grid_idx::CellIdx{N} in one(CellIdx{N}):vsize(ips.grid)
        if ips.grid[grid_idx] in ips.constraint.dest_types
            ips.search_frontier_new_potentials[grid_idx] = UInt32(0)
            enqueue!(ips.search_frontier_buffer, grid_idx)
            LOG_INFERENCE_PATH && println("\tDest at ", grid_idx)
        end
    end

    # Apply the search frontier to the potential field,
    #    each time queueing up neighboring cells to be re-evaluated.
    while !isempty(ips.search_frontier_buffer)
        #TODO: Popping from a Dict is actually super expensive, do something else (e.g. a stack plus this dict, and let the dict grow full)
        grid_idx = dequeue!(ips.search_frontier_buffer)
        new_value = ips.search_frontier_new_potentials[grid_idx]
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

                    current_best_potential = ips.search_frontier_new_potentials[neighbor_grid_idx]
                    already_in_queue::Bool = (current_best_potential != typemax(UInt32))
                    if !already_in_queue
                        current_best_potential = ips.potential[neighbor_grid_idx]
                    end

                    LOG_INFERENCE_PATH && begin
                        println(
                            "\t\tChecking neighbor at ", neighbor_grid_idx, " | value=", ips.grid[neighbor_grid_idx],
                                " | officialPotential=", ips.potential[neighbor_grid_idx],
                                " | upcomingPotential=", already_in_queue ? current_best_potential : "NULL",
                                " | ", ips.grid[neighbor_grid_idx] in ips.constraint.path_types,
                                    "/", neighbor_new_value < current_best_potential
                        )
                    end

                    if (ips.grid[neighbor_grid_idx] in union(ips.constraint.path_types, ips.constraint.source_types)) &&
                        (neighbor_new_value < current_best_potential)
                    #begin
                        if !already_in_queue
                            enqueue!(ips.search_frontier_buffer, neighbor_grid_idx)
                        end
                        ips.search_frontier_new_potentials[neighbor_grid_idx] = neighbor_new_value

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

    function AllInference(paths::Vector{InferPath}, temperature::Real)
        paths_by_type = Dict{UInt8, Vector{Int}}()
        for (i, path) in enumerate(paths)
            for type in path.source_types
                push!(
                    get!(() -> Vector{Int}(), paths_by_type, type),
                    i
                )
            end
        end

        return new(paths, paths_by_type,
                   convert(Float32, temperature))
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
        return new(paths, paths_by_type,
                   temp)
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
    buffer2::Vector{Float32}

    function AllInference_State(constraints::AllInference, grid::CellGrid{N}) where {N}
        return new{N}(
            constraints, grid,
            InferPath_State.(constraints.paths, Ref(grid)),
            Vector{UInt8}(), Vector{Float32}()
        )
    end
end

"
Asks inference how much it likes to apply the given rule to the given location.
The returned value can be compared to other potential applications to decide what to do.
Returns `nothing` if the application should not happen at all.

Make sure to check `inference_exists()` before bothering to execute this.
"
function infer_weight(inference::AllInference_State{N},
                      rule::CellRule, at::CellLine{N},
                      rng::PRNG
                     )::Optional{Float32} where {N}
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

    weight_from_temperature::Float32 = rand(rng, Float32) * inference.source.temperature

    # Accumulate inference weights from path constraints.
    weight_from_path::Float32 = 0
    path_penalty::Float32 = 1
    is_illegal = for_each_cell_in_line(at, true) do idx, cell::CellIdx{N}
        cell_value = inference.grid[cell]
        paths = get(inference.source.paths_by_type, cell_value, nothing)
        if exists(paths)
            # Gather all potentials from paths that originate at this color.
            path_states = (inference.paths[i] for i in paths)
            empty!(inference.buffer2)
            append!(inference.buffer2, (p.potential[cell] for p in path_states))

            # If no path potentials lead to this cell, the move is illegal.
            if none(p -> p >= 0, inference.buffer2)
                return true
            end

            # For each path potential leading to this cell, add it to the final weight.
            # For each path potential not leading to this cell, count it as a penalty.
            for (potential, path_state) in zip(inference.buffer2, path_states)
                if potential < 0
                    path_penalty *= path_state.constraint.penalty
                else
                    weight_from_path += (rand(rng, Float32) * path_state.constraint.temperature) +
                                        (path_state.constraint.invert ? -potential : potential)
                end
            end
        end

        # Nothing illegal about this cell.
        return false
    end
    if is_illegal
        return nothing
    end
    weight_from_path *= path_penalty

    return weight_from_path + weight_from_temperature
end
"
Checks the total potential of a grid cell (mainly for debugging/visualization).
Returns a negative value if the cell's color is off all inference paths.
"
function visualize_weight(inference::AllInference_State{N}, at::CellIdx{N}) where {N}
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