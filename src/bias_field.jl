# How the Field Bias should handle rewrite actions outside its pathing tree.
@bp_enum(BiasFieldOutsideCellMode,
    # Rewrite actions that stray outside the pathing tree are strictly forbidden;
    #    reversing the field does not affect this behavior (they are still forbidden).
    normal,

    # As normal, rewrite actions that go outside the field are forbidden.
    # However this DOES flip when the field flips -- actions INSIDE the pathing tree get forbidden!
    flippable,

    # Actions outside the pathing tree are allowed but have the lowest bias (are done last).
    # When the field is flipped, those outside actions have the highest bias instead.
    soft
)
@bp_enum(BiasFieldComboMode,
    average, min, max,
    first, deviation, diff
)

#= """
Applies a bias in terms of distance from certain cells along certain paths.
Refer to the DSL doc to understand all the bells and whistles.
""" =#
@hybrid struct MarkovBiasFieldImpl <: AbstractMarkovBias
    sources::CellTypeSet
    combo::E_BiasFieldComboMode

    live::Bool
    flipped::Bool

    paths::CellTypeSet # (empty if unused)
    outside_cell_handling::E_BiasFieldOutsideCellMode

    anchors::CellTypeSet # (empty if unused)

    randomness::Float32
    scale::Float32
    penalty_scale::Float32
end

# The immutable version is the standard; immutable is for internal parsing.
const MarkovBiasField = MarkovBiasFieldImpl_Immut
Base.show(io::IO, ::Type{MarkovBiasField}) = print(io, "MarkovBiasField")
const MarkovBiasField_Mutable = MarkovBiasFieldImpl_Mut
Base.convert(::Type{MarkovBiasField}, m::MarkovBiasField_Mutable) = MarkovBiasField(
    getfield.(Ref(m), fieldnames(MarkovBiasField))...
)


const TFieldBiasAnchorBuffers{NGrid} = Tuple{
    Array{Bool, NGrid},
    FloodFillAllocations{Vec{NGrid, Int32}}
}
struct MarkovBiasField_State{NGrid, BHasPathCells, BHasAnchors,
                             TAnchorBuffers<:Optional{TFieldBiasAnchorBuffers{NGrid}}}
    distance_field::Array{UInt32, NGrid} # typemax represents "no path"
    anchor_buffers::TAnchorBuffers
    largest_dist::UInt32
    biases_buffer::Vector{Float32}
end

"""
Returns the largest distance value in the field
  (not including `typemax(UInt32)` which means 'unreachable').

Returns `typemax(UInt32)` if there are no source cells in the entire grid.
"""
function rebuild_distance_field(field::MarkovBiasField,
                                state::MarkovBiasField_State{NGrid, BHasPathCells, BHasAnchors},
                                grid::CellGrid{NGrid},
                                # Convenient type deduction of the integer vector type:
                                ::V = zero(Vec{NGrid, Int32})
                               )::UInt32 where {NGrid, BHasPathCells, BHasAnchors, V<:Vec{NGrid, Int32}}
    grid_idcs_iter = (V(Tuple(i)...) for i in eachindex(CartesianIndices(grid)))

    fn_check_path = if BHasPathCells
        v -> (grid[v] in field.paths)
    else
        v -> true
    end
    fn_check_anchor = if BHasAnchors
        anchor_idcs = Iterators.filter(v -> (grid[v] in field.anchors), grid_idcs_iter)
        function path_connections(v::V, output_list::AbstractVector{v2i})
            function try_side(_axis::Integer, _dir_bool::Integer)
                axis = convert(Int32, _axis)
                dir_bool = convert(Int32, _dir_bool)

                dir_sign::Int32 = Int32.((-1, 1))[dir_bool]
                axis_pos::Int32 = v[axis] + dir_sign
                v2::V = @set v[axis] = axis_pos

                edge::Int32 = Int32.((1, size(grid, _axis)))[dir_bool]
                within_edge::Bool = cmp(axis_pos, edge) != -1

                # Check that the end is a path cell.
                # Ignore the start! We're sometimes coming from a source/anchor cell, not a path one.
                if within_edge && fn_check_path(grid[v2])
                    push!(output_list, v2)
                end
                return nothing
            end
            for axis in Int32(1):Int32(NGrid), s in (0, 1)
                try_side(axis, s)
            end
        end

        (anchor_mask, flood_allocs) = state.anchor_buffers
        fill!(anchor_mask, false)
        flood_fill(
            anchor_idcs,
            path_connections,

            () -> nothing,
            v -> (anchor_mask[v] = true),
            () -> nothing,

            flood_allocs
        )

        v -> anchor_mask[v]
    else
        v -> true
    end

    # Initialize the pathing field.
    @logic_logln "Computing pathing field..."
    fill!(state.distance_field, typemax(UInt32))
    any_sources_exist::Bool = false
    for v in grid_idcs_iter
        if grid[v] in field.sources
            state.distance_field[v] = zero(UInt32)
            @logic_logln "\tSource at " v
            any_sources_exist = true
        end
    end
    if !any_sources_exist
        @logic_logln "No sources were found!"
        return typemax(UInt32)
    end

    # Iterate on the pathing field until it stabilizes.
    max_value::UInt32 = 0
    field_changed::Bool = true
    while field_changed
        field_changed = false

        for v in grid_idcs_iter
            # Compare this cell to the ones ahead of it,
            #    and propagate the field in either direction as applicable.
            cell_a = state.distance_field[v]
            is_invalid_a = (cell_a == typemax(UInt32))

            for a in 1:NGrid
                if v[a] < size(grid, a)
                    v2 = @set v[a] += one(Int32)
                    cell_b = state.distance_field[v2]
                    is_invalid_b = (cell_b == typemax(UInt32))

                    if fn_check_anchor(v2) && fn_check_path(v2)
                        if !is_invalid_a && (is_invalid_b || (cell_a + one(UInt32) < cell_b))
                            @logic_logln "\t\tPush " v.data "(" convert(Int, cell_a) ") to " v2.data
                            field_changed = true
                            cell_b = cell_a + one(UInt32)

                            state.distance_field[v2] = cell_b
                            max_value = max(max_value, cell_b)
                        elseif !is_invalid_b && (is_invalid_a || (cell_b + one(UInt32) < cell_a))
                            @logic_logln "\t\tPush " v2.data "(" convert(Int, cell_b) ") to " v.data
                            field_changed = true
                            cell_a = cell_b + one(UInt32)

                            state.distance_field[v] = cell_a
                            max_value = max(max_value, cell_a)
                        end
                    end
                end
            end
        end
    end
    @logic_logln "\tFinal distance field: " convert.(Ref(Int), state.distance_field) "\n"

    return max_value
end

function markov_bias_state_type(f::MarkovBiasField, ::Type{<:CellGrid{NGrid}},
                                ::PRNG, ::MarkovBiasContext
                               )::Type{<:MarkovBiasField_State} where {NGrid}
    return MarkovBiasField_State{NGrid, !isempty(f.paths), !isempty(f.anchors),
                                 isempty(f.anchors) ? Nothing : TFieldBiasAnchorBuffers{NGrid}}
end
function markov_bias_initialize(f::MarkovBiasField, ::Type{TState},
                                grid::CellGrid{NGrid}, ::PRNG,
                                context::MarkovBiasContext
                               ) where {NGrid, BHasPathCells, BHasAnchors,
                                        TState<:MarkovBiasField_State{NGrid, BHasPathCells, BHasAnchors}}
    # We have a lot to allocate, so dispatch on the allocator type.
    function run_with_alloc(@specialize(alloc))
        V = Vec{NGrid, Int32}
        state = TState(
            markov_allocator_acquire_array(alloc, size(grid), UInt32),
            BHasAnchors ? tuple(
                markov_allocator_acquire_array(alloc, size(grid), Bool),
                FloodFillAllocations(
                    markov_allocator_acquire_set(alloc, V),
                    markov_allocator_acquire_array(alloc, tuple(256), V)
                )
            ) : nothing,
            zero(UInt32),
            markov_allocator_acquire_array(alloc, (), Float32)
        )
        # Setfield package isn't working properly here, for some reason.
        state = TState(
            state.distance_field,
            state.anchor_buffers,
            rebuild_distance_field(f, state, grid),
            state.biases_buffer
        )
        return state
    end
    return run_with_alloc(context.allocator)
end
function markov_bias_cleanup(f::MarkovBiasField, s::MarkovBiasField_State, ctx::MarkovBiasContext)
    # We have a lot to release, so dispatch on the allocator type.
    function run_with_alloc(@specialize(alloc))
        markov_allocator_release_array(alloc, s.distance_field)
        if exists(s.anchor_buffers)
            markov_allocator_release_array(alloc, s.anchor_buffers[1])
            markov_allocator_release_array(alloc, s.anchor_buffers[2].interesting_nodes)
            markov_allocator_release_set(alloc, s.anchor_buffers[2].visited_nodes)
        end
        return nothing
    end
    return run_with_alloc(ctx.allocator)
end

function markov_bias_update(field::MarkovBiasField, state::MarkovBiasField_State{N, BPath, BAnchors},
                            grid::CellGrid{N}, subset::BoxI{N}, old_subset_values::CellGrid{N}
                           )::Nothing where {N, BPath, BAnchors}
    if field.live
        # First check whether the changed area had any relevant cell types.
        involved_cell_types = union(CellTypeSet(old_subset_values),
                                    (grid[v] for v in min_inclusive(subset):max_inclusive(subset)))
        relevant_cell_types = union(field.sources, field.paths, field.anchors)
        if !isempty(intersect(involved_cell_types, relevant_cell_types))
            @set! state.largest_dist = rebuild_distance_field(field, state, grid)
        end
    end
    return state
end
function markov_bias_calculate(field::MarkovBiasField, state::MarkovBiasField_State{N, BPath, BAnchors},
                               grid::CellGrid{N}, at::Union{CellLine{N}, CellRegion{N}},
                               rng::PRNG,
                               # Convenient type-deduction of the integer vector:
                               ::V = zero(Vec{N, Int32})
                              )::Optional{Float32} where {N, BPath, BAnchors, V<:Vec{N, Int32}}
    # First check if this move is legal.
    (on_path, off_path) = if BPath
        let on_path = Ref(false),
            off_path = Ref(false)
          for_each_cell(at) do local_idx, global_idx
              if state.distance_field[global_idx] == typemax(UInt32)
                  off_path[] = true
              else
                  on_path[] = true
              end
          end
          (on_path[], off_path[])
        end
    else
        (true, false)
    end
    # This move is illegal if we're 100% outside the pathing area and outside cells are forbidden.
    if off_path && !on_path &&
        (field.outside_cell_handling == BiasFieldOutsideCellMode.flippable ||
         (field.outside_cell_handling == BiasFieldOutsideCellMode.normal && !field.flipped))
    #begin
        return nothing
    end
    # This move is illegal if we're 100% inside the pathing area and inside cells are forbidden.
    if on_path && !off_path && field.flipped && (field.outside_cell_handling == BiasFieldOutsideCellMode.flippable)
        return nothing
    end
    # Otherwise the move is legal!

    # If randomness is maxed-out, return a uniform-random bias.
    # This skips a lot of work and prevents float issues below.
    if field.randomness >= one(Float32)
        return rand(rng, Float32) * state.largest_dist * field.scale
    end

    # Compute the overall bias.
    function pixel_bias(pos::V)::Float32
        #NOTE: for numeric reasons the bias will approach but never be equal to 0.
        #      Real distance field values will span 1-(max_len+1),
        #        while cells off of the path tree are given fractional values.
        d = state.distance_field[pos]
        return convert(Float32,
            if d == typemax(UInt32)
                # We already know the move is legal; just decide whether to flip its bias.
                bias_flips::Bool = if field.outside_cell_handling in (BiasFieldOutsideCellMode.flippable, BiasFieldOutsideCellMode.soft)
                    true
                elseif field.outside_cell_handling in (BiasFieldOutsideCellMode.normal, )
                    false
                else
                    error("Unhandled: ", field.outside_cell_handling)
                end

                # Now compute the output.
                # Remember that higher values mean higher priority.
                if bias_flips && field.flipped
                    state.largest_dist * field.penalty_scale
                else
                    1.0f0 / field.penalty_scale
                end
            else
                # Remember that higher values mean higher priority.
                if field.flipped
                    d
                else
                    (state.largest_dist - d) + one(UInt32)
                end
            end
        )
    end
    bias::Float32 = if field.combo == BiasFieldComboMode.first
        if at isa CellLine{N}
            pixel_bias(at.start_cell)
        elseif at isa CellRegion{N}
            pixel_bias(min_inclusive(at.b))
        else
            error("Unhandled: ", typeof(at))
        end
    else
        empty!(state.biases_buffer)
        for_each_cell(at) do local_idx, global_idx
            push!(state.biases_buffer, pixel_bias(global_idx))
        end
        if field.combo == BiasFieldComboMode.average
            mean(state.biases_buffer)
        elseif field.combo == BiasFieldComboMode.min
            minimum(state.biases_buffer)
        elseif field.combo == BiasFieldComboMode.max
            maximum(state.biases_buffer)
        elseif field.combo == BiasFieldComboMode.deviation
            # std can at most be half the range, so the value should be doubled for normalization.
            # However in practice it'll be much less, except when some pixels are outside of the path,
            #   so also provide an exponential weighting.
            s = std(state.biases_buffer) * Float32(2)
            if s < Float32(state.largest_dist) # Pixels outside the path tree will have overlarge values
                STD_CURVE = 0.52f0
                s = ((s / Float32(state.largest_dist)) ^ STD_CURVE) * state.largest_dist
            end
            s
        elseif field.combo == BiasFieldComboMode.diff
            # It can already span the whole range, however tends to have a smaller value.
            s = Float32(maximum(state.biases_buffer) - minimum(state.biases_buffer))
            if s < Float32(state.largest_dist) # Pixels outside the path tree will have overlarge values
                DIFF_CURVE = lerp(
                    # Diff starts out like stddev, but can fairly quickly converge on the full range.
                    0.52f0, 1.0f0,
                    clamp(length(state.biases_buffer) / 6.0f0, 0.0f0, 1.0f0)
                )
                s = ((s / Float32(state.largest_dist)) ^ DIFF_CURVE) * state.largest_dist
            end
            s
        else
            error("Unhandled: ", field.combo)
        end
    end

    # Randomize it as necessary.
    if field.randomness > zero(Float32)
        # To get linear weighted randomness,
        #   in this setup where each rewrite move is considered in isolation,
        #   return pow(uniform_f, 1/weight).
        # Bias is integer and can be 0, so here we use the exponent 1/(weight+1) instead.
        weight_curve = 1.0f0 / (bias + 1.0f0)

        # We want a smooth transition from nonrandom (enormously heavy weights, curve approaches +Inf),
        #    to linear weighted-random (curve is unchanged),
        #    to uniform-random (extremely light weights, curve approaches 0.0).
        # We already handled randomness <= 0 and >= 1, so here we can ignore extremes.
        weight_curve = if field.randomness < 0.5f0
            weight_curve / (field.randomness * 2.0f0)
        else
            weight_curve * (1.0f0 - inv_lerp(0.5f0, 1.0f0, field.randomness))
        end

        bias = (rand(rng, Float32) ^ weight_curve)

        # Preserve the magnitude range of 0 - max_path_length.
        bias *= state.largest_dist
    end

    return bias * field.scale
end

function parse_markovjunior_bias(::Val{:field}, inputs::MacroParserInputs,
                                 loc::Optional{LineNumberNode},
                                 expr_args)
    # Construct a default object and customize it with new args as we go.
    field = MarkovBiasField_Mutable(
        CellTypeSet(), BiasFieldComboMode.average,
        false, false,
        CellTypeSet(),
        BiasFieldOutsideCellMode.normal,
        CellTypeSet(),
        0.0f0, 1.0f0, 2.0f0
    )

    # The first argument should be the cell types.
    if isempty(expr_args)
        raise_parse_error(loc, inputs, "Expected at least one argument, for the field's colors!")
    end
    function parse_cell_set(set_expr, msg_name...)::CellTypeSet
        return with_parser_stacktrace(inputs, "Color set \"$(msg_name...)\"") do
            try
                CellTypeSet(set_expr::Symbol)
            catch e
                raise_parse_error(loc, inputs, sprint(showerror, e))
            end
        end
    end
    with_parser_stacktrace(inputs, "Main argument ($(expr_args[1]))") do
        e = expr_args[1]

        eSources = nothing
        ePaths = nothing
        eAnchors = nothing
        if @capture(e, eSources_Symbol->ePaths_Symbol & eAnchors_Symbol) ||
           @capture(e, eSources_Symbol->ePaths_Symbol) ||
           @capture(e, eSources_Symbol)
        #begin
            field.flipped = false
        elseif @capture(e, -eSources_Symbol->ePaths_Symbol & eAnchors_Symbol) ||
               @capture(e,  eSources_Symbol<-ePaths_Symbol & eAnchors_Symbol) ||
               @capture(e, -eSources_Symbol->ePaths_Symbol) ||
               @capture(e,  eSources_Symbol<-ePaths_Symbol) ||
               @capture(e, -eSources_Symbol)
        #begin
            field.flipped = true
        else
            raise_parse_error(loc, inputs, "Expected syntax like 'R', 'R->G', or 'R->G & B'")
        end

        field.sources = parse_cell_set(eSources, "sources")
        exists(ePaths) && (field.paths = parse_cell_set(ePaths, "paths"))
        exists(eAnchors) && (field.anchors = parse_cell_set(eAnchors, "anchors"))
    end
    @markovjunior_assert(!isempty(field.sources))

    # Parse subsequent arguments.
    for arg in @view(expr_args[2:end])
        if arg == :live
            if field.live
                raise_parse_error(loc, inputs, "'live' given more than once")
            else
                field.live = true
            end
        elseif @capture(arg, soft=x_) || (arg == :soft)
            if exists(x) && !isa(x, Real)
                raise_parse_error(loc, inputs, "Value in 'soft=X' should be a number; got ", x)
            elseif isempty(field.paths)
                raise_parse_error(loc, inputs, "'soft' does nothing without path cells")
            elseif field.outside_cell_handling == BiasFieldOutsideCellMode.soft
                raise_parse_error(loc, inputs, "'soft' given more than once")
            elseif field.outside_cell_handling == BiasFieldOutsideCellMode.flippable
                raise_parse_error(loc, inputs, "'soft' and 'forbidden' are mutually incompatible")
            else
                field.outside_cell_handling = BiasFieldOutsideCellMode.soft
                exists(x) && (field.penalty_scale = convert(Float32, x))
            end
        elseif arg == :forbidden
            if isempty(field.paths)
                raise_parse_error(loc, inputs, "'forbidden' does nothing without path cells")
            elseif field.outside_cell_handling == BiasFieldOutsideCellMode.flippable
                raise_parse_error(loc, inputs, "'forbidden' given more than once")
            elseif field.outside_cell_handling == BiasFieldOutsideCellMode.soft
                raise_parse_error(loc, inputs, "'forbidden' and 'soft' are mutually incompatible")
            else
                field.outside_cell_handling = BiasFieldOutsideCellMode.flippable
                field.flipped = true
            end
        elseif arg == :random
            if !iszero(field.randomness)
                raise_parse_error(loc, inputs, "'random'/'randomness=X' given more than once")
            else
                field.randomness = 0.5f0
            end
        elseif @capture(arg, randomness=x_)
            if !iszero(field.randomness)
                raise_parse_error(loc, inputs, "'random'/'randomness=X' given more than once")
            elseif !isa(x, Real)
                raise_parse_error(loc, inputs, "'randomness=X' should be a number from 0 to 1; got ", typeof(x))
            elseif (x < 0.0) || (x > 1.0)
                raise_parse_error(loc, inputs, "'randomness' must be between 0 and 1; got ", x)
            else
                field.randomness = convert(Float32, x)
            end
        elseif @capture(arg, scale=x_)
            if !isone(field.scale)
                raise_parse_error(loc, inputs, "'scale=X' given more than once")
            elseif !isa(x, Real)
                raise_parse_error(loc, inputs, "scale=X' should be a number; got ", typeof(x))
            else
                field.scale = convert(Float32, x)
            end
        elseif @capture(arg, combo=x_)
            if x == :average
                field.combo = BiasFieldComboMode.average
            elseif x == :min
                field.combo = BiasFieldComboMode.min
            elseif x == :max
                field.combo = BiasFieldComboMode.max
            elseif x == :first
                field.combo = BiasFieldComboMode.first
            elseif x == :deviation
                field.combo = BiasFieldComboMode.deviation
            elseif x == :diff
                field.combo = BiasFieldComboMode.diff
            else
                raise_parse_error(loc, inputs, "Unsupported combo! Should be one of ", BiasFieldComboMode.instances())
            end
        else
            raise_parse_error(loc, inputs, "Unsupported argument: \"", arg, "\"")
        end
    end

    return convert(MarkovBiasField, field)
end

dsl_string(f::MarkovBiasField) = string("field(",
    # Source cells:
    (f.flipped && isempty(f.paths)) ? "-" : "",
    dsl_string(f.sources),

    # Path cells:
    isempty(f.paths) ? "" : (
        f.flipped ? "<-" : "->"
    ),
    dsl_string(f.paths),

    # Anchors:
    isempty(f.anchors) ? "" : (
        " & "
    ),
    dsl_string(f.anchors),

    # Now the other arguments.

    # Live:
    f.live ? ", live" : "",

    # Forbidden/soft:
    if f.outside_cell_handling == BiasFieldOutsideCellMode.normal
        ""
    elseif f.outside_cell_handling == BiasFieldOutsideCellMode.flippable
        # Only happens if 'forbidden' was used
        @bp_check(f.flipped)
        ", forbidden"
    elseif f.outside_cell_handling == BiasFieldOutsideCellMode.soft
        if f.penalty_scale == 2.0f0
            ", soft"
        else
            ", soft=$(f.penalty_scale)"
        end
    else
        error("Unhandled: ", f.outside_cell_handling)
    end,

    # Random/randomness:
    if f.randomness == 0.5f0
        ", random"
    elseif !iszero(f.randomness)
        ", randomness=$(f.randomness)"
    else
        ""
    end,

    # Scale:
    if !isone(f.scale)
        ", scale=$(f.scale)"
    else
        ""
    end,

    # Combo (for simplicity, always specify even if default):
    ", combo=", f.combo,

    ")"
)