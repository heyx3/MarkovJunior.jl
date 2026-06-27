####################################
#  Rewrite rule: individual Cells

const RewriteRuleCell_Set = CellTypeSet
const RewriteRuleCell_List = UpTo{N_CELL_TYPES, UInt8}
struct RewriteRuleCell_Wildcard end
struct RewriteRuleCell_Lookup{TIndex}
    source_idx::TIndex
end

"The source entries for a rewrite rule can be a single value, ordered set, or a wildcard"
const RewriteRuleCellSource = Union{UInt8, RewriteRuleCell_Set, RewriteRuleCell_Wildcard}
"The destination entries for a rewrite rule can be a single value, unordered set, list (matching source set size), or a wildcard"
const RewriteRuleCellDest{TIndex} = Union{UInt8, RewriteRuleCell_Set,
                                          RewriteRuleCell_List, RewriteRuleCell_Wildcard,
                                          RewriteRuleCell_Lookup{TIndex}}

const RewriteCell{TIndex} = Tuple{RewriteRuleCellSource, RewriteRuleCellDest{TIndex}}

match_rewrite_source(rule::UInt8, value::UInt8) = (rule == value)
match_rewrite_source(rule::RewriteRuleCell_Set, value::UInt8) = (value in rule)
match_rewrite_source(::RewriteRuleCell_Wildcard, ::UInt8) = true

# 1D lookups provide the source values as a Tuple.
pick_rewrite_value(dest::UInt8,                       src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = dest
pick_rewrite_value(dest::RewriteRuleCell_Set,         src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = rand(rng, dest)
pick_rewrite_value(dest::RewriteRuleCell_List,        src::RewriteRuleCell_Set,   src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = dest[cell_set_index_of(src, src_values[self_idx])]
pick_rewrite_value(dest::RewriteRuleCell_Lookup{Int}, src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = src_values[dest.source_idx]
pick_rewrite_value(dest::RewriteRuleCell_Wildcard,    src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = src_values[self_idx]

# Multi-dimensional lookups provide the source values as some array-view
#    (correcting for symmetry, meaning the first array axis is the first rule axis and so on).
pick_rewrite_value(dest::UInt8,                                  src::RewriteRuleCellSource, src_values::AbstractArray{UInt8, N}, self_idx::NTuple{N, Int}, rng::PRNG) where {N} = dest
pick_rewrite_value(dest::RewriteRuleCell_Set,                    src::RewriteRuleCellSource, src_values::AbstractArray{UInt8, N}, self_idx::NTuple{N, Int}, rng::PRNG) where {N} = rand(rng, dest)
pick_rewrite_value(dest::RewriteRuleCell_List,                   src::RewriteRuleCell_Set,   src_values::AbstractArray{UInt8, N}, self_idx::NTuple{N, Int}, rng::PRNG) where {N} = dest[cell_set_index_of(src, src_values[self_idx...])]
pick_rewrite_value(dest::RewriteRuleCell_Lookup{NTuple{N, Int}}, src::RewriteRuleCellSource, src_values::AbstractArray{UInt8, N}, self_idx::NTuple{N, Int}, rng::PRNG) where {N} = src_values[dest.source_idx...]
pick_rewrite_value(dest::RewriteRuleCell_Wildcard,               src::RewriteRuleCellSource, src_values::AbstractArray{UInt8, N}, self_idx::NTuple{N, Int}, rng::PRNG) where {N} = src_values[self_idx...]


###########################
#  Rewrite rule: 1D strip

"A description of extra allowed symmetries starting at a specific axis, optionally different starting axes per-direction"
const RewriteRule_TailSymmetry = Union{Nothing, Int, NTuple{2, Optional{Int}}}

struct RewriteRule_Strip{NCells, TCells<:NTuple{NCells, RewriteCell{Int}}}
    cells::TCells
    mask::Union{Nothing, Float32, NTuple{2, Float32}} # Disabled, constant, or random range
    weight::Float32
    explicit_symmetries::Vector{GridDir}
    tail_symmetry::RewriteRule_TailSymmetry
end
Base.:(==)(a::RewriteRule_Strip{N, T}, b::RewriteRule_Strip{N, T}) where {N, T} = (
    a.cells == b.cells &&
    a.mask == b.mask &&
    a.weight == b.weight &&
    a.explicit_symmetries == b.explicit_symmetries &&
    a.tail_symmetry == b.tail_symmetry
)

const MaskGrid{N} = Array{Float32, N}

pick_mask(r::RewriteRule_Strip, rng::PRNG)::Float32 = pick_mask_impl(r.mask, rng)
pick_mask_impl(::Nothing, ::PRNG) = 1.0f0
pick_mask_impl(f::Float32, ::PRNG) = f
pick_mask_impl((a, b)::NTuple{2, Float32}, rng::PRNG) = lerp(a, b, rand(rng, Float32))

"
Checks for a rule matching against the given grid when applied to the given cell.
Assumes for performance that the rule does fit into the grid.
"
function rule_matches(r::RewriteRule_Strip{NCells, TCells}, grid::CellGrid{NDims},
                      c_start::CellIdx{NDims}, dir::GridDir
                     )::Bool where {NCells, TCells, NDims}
    @markovjunior_assert begin
        c_end = grid_dir_pos_along(dir, c_start, NCells - 1)
        grid_range = one(VecI{NDims}):vsize(grid)
        (c_start in grid_range) && (c_end in grid_range)
    end "Rule is outside grid bounds! $c_start along $dir by $(NCells-1) vs $(vsize(grid))"

    return all(ntuple(Val(NCells)) do i
        c_next = grid_dir_pos_along(dir, c_start, i-1)
        return match_rewrite_source(r.cells[i][1], grid[c_next])
    end)
end

"
Finds every match for a given rule that touches a given subset of the grid,
  excluding rule placements that start at cells which fail the given mask.

For each potential match, your lambda receives `(start_cell, grid_dir, match_succeeded)`.
Note this does not include cells which are outside the mask.

If your lambda returns anything other than `nothing`,
  the function immediately stops checking for matches and returns that value;
  otherwise we return `nothing` at the end of the search.
"
function visit_rule_match_data(process_candidate::TLambda,
                               r::RewriteRule_Strip{NCells, TCells}, grid::CellGrid{NDims},
                               grid_mask::Optional{MaskGrid{NDims}}, rule_chosen_mask::Float32,
                               inner_subset::Bplus.Math.BoxI{NDims}
                              )::Optional where {NCells, TCells, NDims, TLambda}
    function process_dir(dir::GridDir)
        clamp_range = (
            if dir.sign > 0
                1
            else
                NCells
            end,
            size(grid, dir.axis) + if dir.sign > 0
                -(NCells - 1)
            else
                0
            end
        )

        first_pos = min_inclusive(inner_subset)
        if dir.sign > 0
            first_pos = grid_dir_pos_along(dir, first_pos, -(NCells - 1))
        end
        @set! first_pos[dir.axis] = clamp(first_pos[dir.axis], clamp_range...)

        last_pos = max_inclusive(inner_subset)
        if dir.sign < 0
            last_pos = grid_dir_pos_along(dir, last_pos, -(NCells - 1))
        end
        @set! last_pos[dir.axis] = clamp(last_pos[dir.axis], clamp_range...)

        # Note that if first_pos isn't completely behind last_pos, it means
        #    the rule can't fit in the grid along this direction.
        for rule_start_pos::CellIdx{NDims} in first_pos:last_pos
            if isnothing(grid_mask) || grid_mask[rule_start_pos] <= rule_chosen_mask
                user_out = process_candidate(rule_start_pos, dir,
                                             rule_matches(r, grid, rule_start_pos, dir))
                exists(user_out) && return user_out
            end
        end
    end

    for dir in r.explicit_symmetries
        uo = process_dir(dir)
        exists(uo) && return uo
    end

    if isnothing(r.tail_symmetry)
        # No extra symmetries to process.
    elseif r.tail_symmetry isa Int
        for axis in r.tail_symmetry:NDims
            for sign in (-1, 1)
                uo = process_dir(GridDir(axis, sign))
                exists(uo) && return uo
            end
        end
    elseif r.tail_symmetry isa NTuple{2, Optional{Int}}
        (a_neg, a_pos) = r.tail_symmetry
        exists(a_neg) && for axis in a_neg:NDims
            uo = process_dir(GridDir(axis, -1))
            exists(uo) && return uo
        end
        exists(a_pos) && for axis in a_pos:NDims
            uo = process_dir(GridDir(axis, +1))
            exist(uo) && return uo
        end
    else
        error("Unhandled: ", typeof(r.tail_symmetry))
    end

    return nothing
end


####################################
#  Rewrite rule: multidimensional

# Note that some definitions from the 1D "strip" rewrite rules are also used here.

"
Data used to generate the symmetries of a multi-dimensional rewrite rule.
Multi-dimensional rule orientation is conceived as a list of choices --
  for each rule axis, pick a corresponding grid axis (positive and/or negative).
"
struct RewriteRule_MD_Symmetry_Definition
    # Each entry is a pair:
    #   1. Possibile orientations for a small group of rule axes, as a Matrix
    #   2. If only one rule axis is involved, this is its tail symmetry. Otherwise this is `nothing`.
    # The first matrix column tells you the rule axes involved -- each row corresponds to that axis.
    # Subsequent columns represent allowed orientations of rule axes along the grid.
    # Each axis orientation is a signed integer --
    #    +a means "along positive axis a"; -a means "along negative axis a".
    # The tail symmetry for each axis tells you whether arbitrary extra dimensions are allowed to be chosen.
    # Any rule axes not mentioned here are allowed to take on any orientation.
    grid_axis_choices::Vector{Pair{Matrix{Int}, RewriteRule_TailSymmetry}}

    # Each set of rule axes listed here may not flip relative to each other,
    #   i.e. they may only be rotated into their grid orientation.
    chiral_groups::Vector{Set{Int}}
end
Base.:(==)(a::RewriteRule_MD_Symmetry_Definition, b::RewriteRule_MD_Symmetry_Definition) = (
    a.grid_axis_choices == b.grid_axis_choices &&
    a.chiral_groups == b.chiral_groups
)
"The empty-constructor creates the maximum possible symmetry"
RewriteRule_MD_Symmetry_Definition() = RewriteRule_MD_Symmetry_Definition(
    Vector{Pair{Matrix{Int}, RewriteRule_TailSymmetry}}(),
    Vector{Set{Int}}()
)

log_md_symmetry_logic() = false
macro md_symm_log(args...)
    MJ = @__MODULE__
    return :(
        if $MJ.log_md_symmetry_logic()
            print(stderr, $(esc.(args)...))
        end
    )
end
macro md_symm_logln(args...)
    MJ = @__MODULE__
    return :(
        if $MJ.log_md_symmetry_logic()
            println(stderr, $(esc.(args)...))
        end
    )
end
macro md_symm_display(args...)
    MJ = @__MODULE__
    return :(
        if $MJ.log_md_symmetry_logic()
            display($(esc.(args)...))
        end
    )
end

"
Generates all legal orientations of the given rule block in the given grid.
Returns a matrix where each column is an orientation -- mapping each rule axis (row) to a grid axis (element).
"
function find_all_md_symmetries(def::RewriteRule_MD_Symmetry_Definition,
                                n_rule_dims::Int, n_grid_dims::Int
                               )::AbstractMatrix{GridDir}
    @md_symm_logln("Finding symmetries for definition with ",
                    n_rule_dims, " Rule dims and ", n_grid_dims, " Grid dims")
    md_symm_loggable_dirs(dirs) = Iterators.flatten(Iterators.map(dirs) do d::GridDir
        return (d.axis == -1) ? "__" : d.sign * d.axis
    end)

    options = Vector{Vector{GridDir}}() #TODO: Use an ElasticArrays.jl matrix instead
    UNCHOSEN_GRID_DIR = GridDir(-1, -1, true)
    FRESH_OPTION = map(i->UNCHOSEN_GRID_DIR, 1:n_rule_dims)

    # Occasionally we track invalid options and then remove them in a second pass.
    options_to_remove = Vector{Union{Int, UnitRange{Int}, StepRange{Int, Int}}}()
    function execute_removals()
        # Iterate from end to start, for performance AND simplicity.
        for bad_i in Iterators.reverse(options_to_remove)
            deleteat!(options, bad_i)
        end
        empty!(options_to_remove)
    end

    set_rule_axes = fill(false, n_rule_dims)
    unused_grid_axes = Vector{Bool}(undef, n_grid_dims) # Buffer used below

    #TODO: Track chiral groups and filter out options which violate them.

    # For each new choice to make, try pairing it against every existing choice.
    first_run::Bool = true
    for (raw_new_orientations, tail_symmetry) in def.grid_axis_choices
        orientation_rule_axes::AbstractVector{Int} = @view raw_new_orientations[:, 1]
        @markovjunior_assert(all(a -> (a>0) && (a<=n_rule_dims), orientation_rule_axes),
                             "Some named rule axes are outside the range (1, ", n_rule_dims, "): ",
                               orientation_rule_axes)
        n_orientation_axes::Int = length(orientation_rule_axes)
        @markovjunior_assert(none(a -> set_rule_axes[a], orientation_rule_axes),
                             "Some rule axes assigned to grid axis in more than one way! At least one of ",
                               orientation_rule_axes)
        for a in orientation_rule_axes
            set_rule_axes[a] = true
        end

        new_orientations::AbstractMatrix{Int} = @view raw_new_orientations[:, 2:end]
        n_new_orientations::Int = size(new_orientations, 2)
        new_orientation_columns::@IterOf{AbstractVector{Int}} = Iterators.map(
            i -> @view(new_orientations[:, i]),
            1:n_new_orientations
        )
        # Incorporate tail symmetry.
        tail_symmetry_columns = Vector{Vector{Int}}() #TODO: Each column is 1-long, so we should use a Vector{Int} and then reshape to matrix
        if exists(tail_symmetry)
            @bp_check(n_orientation_axes == 1,
                      "Can't apply tail symmetry to an axis group! Only to an individual axis")

            if tail_symmetry isa Int
                for tail_grid_axis in tail_symmetry:n_grid_dims
                    for tail_sign in (-1, 1)
                        push!(tail_symmetry_columns, [ tail_grid_axis * tail_sign ])
                    end
                end
            elseif tail_symmetry isa NTuple{2, Optional{Int}}
                exists(tail_symmetry[1]) && for tail_grid_axis in tail_symmetry[1]:n_grid_dims
                    push!(tail_symmetry_columns, [ -tail_grid_axis ])
                end
                exists(tail_symmetry[2]) && for tail_grid_axis in tail_symmetry[2]:n_grid_dims
                    push!(tail_symmetry_columns, [ tail_grid_axis ])
                end
            else
                error("Unhandled: ", typeof(tail_symmetry))
            end

        end
        # *Always* append the tail symmetry, even if it's empty and inapplicable,
        #    to improve type-stability.
        n_new_orientations += length(tail_symmetry_columns)
        new_orientation_columns = Iterators.flatten((
            new_orientation_columns,
            tail_symmetry_columns
        ))
        @md_symm_logln("Next Group: [", iter_join(orientation_rule_axes, ",")..., "]")
        @md_symm_logln("Entries: [")
        for entry in new_orientation_columns
            @md_symm_logln("\t[", iter_join(entry, ", ")..., "]")
        end
        @md_symm_logln("]")

        # If no options are available, then it's either the first iteration or there are no valid symmetries.
        if isempty(options)
            if first_run
                first_run = false
                for new_orientation in new_orientation_columns
                    push!(options, map(identity, FRESH_OPTION))
                    new_option = options[end]

                    for (a_rule::Int, a_orientation::Int) in zip(orientation_rule_axes, new_orientation)
                        @markovjunior_assert(new_option[a_rule].axis == -1,
                                             "Rule axis set more than once: ", a_rule)
                        options[end][a_rule] = GridDir(abs(a_orientation), sign(a_orientation))
                    end

                    #TODO: Check chirality
                    @md_symm_logln("New option: [ ", iter_join(md_symm_loggable_dirs(new_option), ", ")..., " ]")
                end
            else
                break
            end
        # Otherwise, try pairing every existing option with each of these new choices.
        else
            # Copy the original options once for each new choice.
            n_original_options = length(options)
            resize!(options, n_original_options * n_new_orientations)
            for src_option_i in 1:n_original_options
                for permutation_i in 2:n_new_orientations # Start at 2 to leave the original options alone
                    dest_option_i = src_option_i + (n_original_options * (permutation_i - 1))
                    options[dest_option_i] = copy(options[src_option_i])
                end
            end

            # Now apply each permutation to all its original options.
            for (permutation_i, new_orientation) in enumerate(new_orientation_columns)
                for src_option_i in 1:n_original_options
                    dest_option_i = src_option_i + (n_original_options * (permutation_i - 1))
                    option = options[dest_option_i]
                    @md_symm_logln(
                        "Applying permutation [",
                        iter_join(new_orientation, ",")...,
                        "] to option [",
                            iter_join(md_symm_loggable_dirs(option), ",")...,
                        "]..."
                    )

                    # Apply this permutation of choices to the rule axes.
                    for (a_rule::Int, a_orientation::Int) in zip(orientation_rule_axes, new_orientation)
                        @markovjunior_assert(option[a_rule].axis == -1,
                                             "More than one symmetry statement is writing to rule axis ", a_rule)
                        # If any other rule axes aleady map to this same grid axis, this option is not viable.
                        if any(a -> option[a].axis == abs(a_orientation), 1:n_rule_dims)
                            @md_symm_logln(
                                "\tThis option is invalid! Rule axis ", a_rule,
                                " will map to grid axis ", abs(a_orientation),
                                ", but it was already being used"
                            )
                            push!(options_to_remove, dest_option_i)
                            break
                        end

                        option[a_rule] = GridDir(abs(a_orientation), sign(a_orientation))
                    end
                    if !isempty(options_to_remove) && (options_to_remove[end] == dest_option_i)
                        continue
                    end
                    @md_symm_logln(
                        "\tNow the option is [",
                        iter_join(md_symm_loggable_dirs(option), ",")..., "]"
                    )

                    #TODO: Check chirality
                end
            end

            execute_removals()
        end
    end

    # Any unmentioned rule axes are implicitly unconstrained -- they can orient along any world axis.
    for a_rule in 1:n_rule_dims if !set_rule_axes[a_rule]
        @md_symm_logln(
            "Rule axis ", a_rule,
            " was not mentioned anywhere so it has implicit tail symmetry of `1...`."
        )

        # If there are no explicit symmetries/existing options at all,
        #    then this rule axis provides an initial set of options.
        if count(set_rule_axes) < 1
            for a_grid in 1:n_grid_dims
                for s_grid in (-1, +1)
                    option = map(i -> UNCHOSEN_GRID_DIR, 1:n_rule_dims)
                    option[a_rule] = GridDir(a_grid, s_grid)
                    push!(options, option)
                end
            end
        else
            # The specific grid axes available will be different for each option,
            #    but the number of choices should be constant.
            n_set_rule_axes = count(set_rule_axes)
            n_original_options = length(options)
            n_grid_axis_choices = n_grid_dims - n_set_rule_axes
            push!(options_to_remove, 1:n_original_options) # After the loop, the original options are replaced
            for src_option_i in 1:n_original_options
                src_option = options[src_option_i]
                @markovjunior_assert(src_option[a_rule].axis == -1,
                                    "Rule-axis ", a_rule, " actually was set??")
                # Gather the candidate grid axes.
                fill!(unused_grid_axes, true)
                for dir in src_option
                    (dir.axis != -1) && (unused_grid_axes[dir.axis] = false)
                end
                @markovjunior_assert(count(unused_grid_axes) == n_grid_axis_choices,
                                    "Expected ", n_grid_axis_choices,
                                    " unused grid axes for each option at this point, but found ",
                                    unused_grid_axes, " in option ", src_option_i)
                # Try using each one.
                for a_grid in 1:n_grid_dims if unused_grid_axes[a_grid]
                    for a_grid_sign in (-1, 1)
                        a_grid_dir = GridDir(a_grid, a_grid_sign)
                        #TODO: Check chirality against this choice

                        option = map(identity, src_option)
                        option[a_rule] = a_grid_dir
                        push!(options, option)
                    end
                end end
            end
        end

        # Finish up.
        execute_removals()
        set_rule_axes[a_rule] = true
    end end

    # If there are no options then the output type can't be deduced automatically :(
    result = isempty(options) ? Matrix{GridDir}(undef, n_rule_dims, 0) : hcat(options...)

    @md_symm_logln("Final symmetry output: [")
    for i_rule_axis in 1:size(result, 1)
        for i_entry in 1:size(result, 2)
            @md_symm_log(md_symm_loggable_dirs((result[i_rule_axis, i_entry], ))..., " ")
        end
        @md_symm_logln()
    end
    @md_symm_logln("]\n")
    return result
end


const RewriteCell_MD{NDims} = RewriteCell{NTuple{NDims, Int}}

"
An optimized lookup of the symmetries for a multidimensional rule embedded in a specific grid.

Iterate over individual orientations from this set with `orientations(o)`.
Release this struct's allocations with `close(o)`.
"
struct RewriteRule_MD_Orientations{NGridDims}
    # Number of availble orientations (i.e. number of columns in the below tables)
    n_entries::Int

    # Each column is an entry; each row is a rule axis; each element maps one rule axis to the grid.
    rule_to_grid::Matrix{GridDir}
    # Each column represents one column from 'rule_to_grid' (i.e. one orientation entry);
    #    each row is one grid axis, storing its optional mapping to the rule axes.
    # If a value's axis is -1, no rule axis is oriented along here.
    # A value's sign means the same thing as the other table --
    #    the direction of rule along this grid axis.
    grid_to_rule::Matrix{GridDir}

    # For each orientation ('n_entries'), this contains the original rule block
    #    transformed by that permutation and with dimensions matching that of the grid.
    # Its contents are also modified accordingly (e.g. "Lookup" cells' indices are permuted).
    rule_permutations::Vector{Array{RewriteCell_MD{NGridDims}, NGridDims}}
    # For each rule permutation, this is an array of identical size,
    #    meant to be a buffer when applying rewrite rules.
    rule_permutation_cells_buffer::Vector{CellGrid{NGridDims}}

    allocator::AbstractMarkovAllocator
end
function RewriteRule_MD_Orientations(rule_array::Array{RewriteCell_MD{NRuleDims}, NRuleDims},
                                     symmetries::RewriteRule_MD_Symmetry_Definition,
                                     ::Val{NGridDims},
                                     allocator::AbstractMarkovAllocator
                                    ) where {NRuleDims, NGridDims}
    @bp_check(NGridDims > 0, "Grid must be at least 1D; got ", NGridDims)
    @bp_check(NRuleDims in 1:NGridDims,
              "Invalid rule dimension (", NRuleDims,
                ") against grid dimension ", NGridDims)

    # Build the initial orientation matrix.
    raw_rule_to_grid = find_all_md_symmetries(symmetries, NRuleDims, NGridDims)
    n_entries::Int = size(raw_rule_to_grid, 2)
    rule_to_grid::Matrix{GridDir} = markov_allocator_acquire_array(allocator, (NRuleDims, n_entries), GridDir)
    copy!(rule_to_grid, raw_rule_to_grid)

    # Build the reverse lookup.
    grid_to_rule::Matrix{GridDir} = markov_allocator_acquire_array(allocator, (NGridDims, n_entries), GridDir)
    fill!(grid_to_rule, GridDir(-1, -1, true))
    for entry_i in 1:n_entries
        for rule_axis in 1:NRuleDims
            dir = rule_to_grid[rule_axis, entry_i]
            grid_to_rule[dir.axis, entry_i] = GridDir(rule_axis, dir.sign)
        end
    end

    # Build the permuted rule matrices.
    PermArrayType = Array{RewriteCell_MD{NGridDims}, NGridDims}
    rule_permutations = markov_allocator_acquire_array(allocator, (n_entries, ), PermArrayType)
    rule_permutation_apply_buffers = markov_allocator_acquire_array(
        allocator,
        (n_entries, ), CellGrid{NGridDims}
    )
    @md_symm_logln("=-=-=-=-=-=-=-\nCrunching orientations for rule:\n", sprint(dump, rule_array),
                   "\n\trule_to_grid: ", sprint(dump, rule_to_grid),
                   "\n\tgrid_to_rule: ", sprint(dump, grid_to_rule))
    for permutation_i in 1:n_entries
        perm_rule_to_grid::AbstractVector{GridDir} = @view rule_to_grid[:, permutation_i]
        perm_grid_to_rule::AbstractVector{GridDir} = @view grid_to_rule[:, permutation_i]
        # Append stub dimensions until the rule matches the world grid.
        rule_full_d = reshape(rule_array, (
            size(rule_array)...,
            (1 for i in (ndims(rule_array)+1):NGridDims)...
        ))
        # Now shuffle the rule axes to follow the permutation.
        #= Sample case:
            rule is [ a, b, c ]
            perm is (-3, )
            rule_to_grid is [ -3 ]
            grid_to_rule is [ 0 0 -1 0 ]
            rule_full_d is [ a, b, c ;;;; ]
            rule_full_permutation should be (2, 3, 1, 4)
        =#
        prev_singleton_rule_axis = Ref(ndims(rule_array))
        axis_permutation::NTuple{NGridDims, Int} = ntuple(Val(ndims(rule_full_d))) do i
            if perm_grid_to_rule[i].axis > 0
                perm_grid_to_rule[i].axis
            else
                prev_singleton_rule_axis[] += 1
                prev_singleton_rule_axis[]
            end
        end
        @md_symm_logln("\t", permutation_i, ": ", ntuple(Val(NRuleDims)) do ar
            perm_rule_to_grid[ar].axis * perm_rule_to_grid[ar].sign
        end, " using permutedims ", axis_permutation, " (from grid side: ", ntuple(Val(NGridDims)) do ag
            if perm_grid_to_rule[ag].axis < 1
                0
            else
                perm_grid_to_rule[ag].axis * perm_grid_to_rule[ag].sign
            end
        end, ")")
        permuted_view = PermutedDimsArray(rule_full_d, axis_permutation)
        # Note that there is no lazy "flip axis" operator, so we wait to do that
        #   until after the permutation is copied into a new array.

        # Finally, generate a new non-lazy array to hold this permutation.
        permuted_allocation = markov_allocator_acquire_array(
            allocator,
            size(permuted_view), RewriteCell_MD{NGridDims}
        )
        rule_permutations[permutation_i] = permuted_allocation
        # Copy in the permuted rule.
        # Certain kinds of cells (i.e. RewriteRuleCell_Lookup) require changes to match the permutation.
        for inds::CartesianIndex{NGridDims} in eachindex(permuted_view)
            cell::RewriteCell_MD{NRuleDims} = permuted_view[inds]
            (cell_src, cell_dest) = cell

            # Process the cell contents.
            cell_dest = if cell_dest isa RewriteRuleCell_Lookup
                RewriteRuleCell_Lookup(ntuple(Val(NGridDims)) do a_grid::Int
                    if perm_grid_to_rule[a_grid].axis < 1
                        # This isn't part of the rule orientation, so it's a stub dimension.
                        1
                    else
                        a_rule            = perm_grid_to_rule[a_grid].axis
                        a_grid_is_flipped = perm_grid_to_rule[a_grid].sign < 0

                        pos_along_axis = cell_dest.source_idx[a_rule]
                        max_along_axis = size(permuted_allocation, a_grid)
                        if a_grid_is_flipped
                            pos_along_axis = max_along_axis - pos_along_axis + 1
                        end

                        pos_along_axis
                    end
                end)
            else
                cell_dest
            end

            permuted_allocation[inds] = (cell_src, cell_dest)
        end
        # Apply the flipping of some axes.
        for i_rule in 1:NRuleDims
            if perm_rule_to_grid[i_rule].sign < 0
                @md_symm_logln("\t\t(flip along rule axis ", i_rule, ")")
                reverse!(permuted_allocation, dims=perm_rule_to_grid[i_rule].axis)
            end
        end
        @md_symm_logln("\t\t", sprint(dump, permuted_allocation))

        # Allocate a parallel array for rule application.
        rule_permutation_apply_buffers[permutation_i] = markov_allocator_acquire_array(
            allocator,
            size(permuted_allocation),
            UInt8
        )
    end

    return RewriteRule_MD_Orientations{NGridDims}(
        n_entries,
        rule_to_grid, grid_to_rule,
        rule_permutations, rule_permutation_apply_buffers,
        allocator
    )
end
Base.close(o::RewriteRule_MD_Orientations) = markov_allocator_release_array.(Ref(o.allocator), (
    o.rule_to_grid,
    o.grid_to_rule,
    o.rule_permutations...,
    o.rule_permutations,
    o.rule_permutation_cells_buffer...,
    o.rule_permutation_cells_buffer
))

"The concrete type of a view into one orientation of a multidimensional rewrite rule"
const RewriteRule_MD_OrientationView = let m = Matrix{GridDir}(undef, 3, 3)
    typeof(@view(m[:, 1]))
end

"A single orientation taken from `RewriteRule_MD_Orientations` (by iterating over `orientations(os)`)"
struct RewriteRule_MD_Orientation{NGridDims}
    # Each field is a column from the corresponding table in the Orientations struct above.
    rule_to_grid::RewriteRule_MD_OrientationView
    grid_to_rule::RewriteRule_MD_OrientationView
    rule_permutation::Array{RewriteCell_MD{NGridDims}, NGridDims}
    application_buffer::Array{UInt8, NGridDims}
end
orientations(os::RewriteRule_MD_Orientations{N}) where {N} = Iterators.map(1:os.n_entries) do i
    return RewriteRule_MD_Orientation{N}(
        @view(os.rule_to_grid[:, i]),
        @view(os.grid_to_rule[:, i]),
        os.rule_permutations[i],
        os.rule_permutation_cells_buffer[i]
    )
end

struct RewriteRule_MD{NDims}
    cells::Array{RewriteCell_MD{NDims}, NDims}
    mask::Union{Nothing, Float32, NTuple{2, Float32}} # Disabled, constant, or random range
    weight::Float32
    symmetry::RewriteRule_MD_Symmetry_Definition
end
Base.:(==)(a::RewriteRule_MD{NDims1}, b::RewriteRule_MD{NDims2}) where {NDims1, NDims2} = (
    NDims1 == NDims2 &&
    a.cells == b.cells &&
    a.mask == b.mask &&
    a.weight == b.weight &&
    a.symmetry == b.symmetry
)
Base.ndims(r::RewriteRule_MD{NDims}) where {NDims} = NDims

pick_mask(r::RewriteRule_MD, rng::PRNG)::Float32 = pick_mask_impl(r.mask, rng)

"
Checks for a rule matching against the given grid when applied starting at the given grid cell.

Assumes for performance that the rule doesn't stick outside the grid.
Additionally, symmetry permutation should already be baked into the rule array
  (e.g. it must have as many dimensions as the grid).
"
function rule_matches(rule_array::Array{RewriteCell_MD{NDims}, NDims},
                      grid_array::CellGrid{NDims},
                      rule_start_in_grid::CellIdx{NDims}
                     )::Bool where {NDims}
    # Assert the rule is in-bounds.
    @markovjunior_assert(
        all((rule_start_in_grid > zero(Int32)) &
            ((rule_start_in_grid + vsize(rule_array) - one(Int32)) <= vsize(grid_array))),
        "Rule starts at ", rule_start_in_grid, " and extends by ", vsize(rule_array),
          ", while grid is ", vsize(grid_array),
          ". Per-axis: min=", rule_start_in_grid > zero(Int32),
                    "  max=", (rule_start_in_grid + vsize(rule_array) - one(Int32)) <= vsize(grid_array)
    )

    # Try the rule.
    @md_symm_logln "Matching a " vsize(rule_array) " rule to the grid at " rule_start_in_grid "..."
    result = all(one(CellIdx{NDims}):vsize(rule_array)) do r_pos
        g_pos = rule_start_in_grid + r_pos - one(Int32)
        cell_rule = rule_array[r_pos][1]
        value = grid_array[g_pos]
        result = match_rewrite_source(cell_rule, value)
        @md_symm_log(
            "    Cell from ", r_pos, " to ", g_pos,
            ", ", cell_rule, " against ", value, ": ", result,
            (r_pos == vsize(rule_array)) ? "\n" : ""
        )
        return result
    end
    @md_symm_logln
    return result
end

"
Finds every match for a given rule that touches a given subset of the grid,
  excluding rule placements that start at cells which fail the given mask.

For each potential match, your lambda receives `(start_cell, orientation, match_succeeded)`.
Note this does not include cells which are outside the mask.

If your lambda returns anything other than `nothing`,
  the function immediately stops checking for matches and returns that value;
  otherwise we return `nothing` at the end of the search.
"
function visit_rule_match_data(process_candidate::TLambda,
                               rule::RewriteRule_MD{NRuleDims},
                               cached_orientations::RewriteRule_MD_Orientations{NGridDims},
                               grid::CellGrid{NGridDims},
                               grid_mask::Optional{MaskGrid{NGridDims}}, rule_chosen_mask::Float32,
                               grid_inner_subset::BoxI{NGridDims}
                              )::Optional where {NRuleDims, NGridDims, TLambda}
    IVec = VecI{NGridDims}
    @md_symm_logln("\n=====-----=======---======\nMatching this rule:\n", sprint(dump, rule.cells),
                   "\n against this grid:\n", sprint(dump, grid))
    for orientation::RewriteRule_MD_Orientation in orientations(cached_orientations)
        perm_size = vsize(orientation.rule_permutation)
        @md_symm_logln("\n\t=======================\nFinding matches for orientation ", Tuple(orientation.rule_to_grid), "...")
        @md_symm_logln("\t\tRule: ", sprint(dump, orientation.rule_permutation))
        start_corner_min::IVec = max(one(IVec),
                                     min_inclusive(grid_inner_subset) - perm_size + one(Int32))
        start_corner_max::IVec = min(vsize(grid) - perm_size + 1,
                                     max_inclusive(grid_inner_subset))
        for start_pos in start_corner_min:start_corner_max
            if isnothing(grid_mask) || (grid_mask[start_pos] <= rule_chosen_mask)
                matches = rule_matches(orientation.rule_permutation, grid, start_pos)
                matches && @md_symm_logln("\t\t\tMATCH at ", start_pos)
                user_out = process_candidate(start_pos, orientation, matches)
                exists(user_out) && return user_out
            end
        end
    end

    return nothing
end


####################################
#  Rewrite rule: cache of matches

const CachedRuleApplication{NDims} = Tuple{
    CellIdx{NDims},
    Union{GridDir, RewriteRule_MD_Orientation{NDims}}
}
"Efficiently tracks all possible applications of rewrite rules to a grid"
struct RewriteCache{NDims, NRules,
                    TGrid<:CellGrid{NDims},
                    TRules<:NTuple{NRules, Union{RewriteRule_Strip, RewriteRule_MD}}}
    grid::TGrid
    rules::TRules
    rule_md_orientations::NTuple{NRules, Optional{RewriteRule_MD_Orientations{NDims}}}

    mask_grid::Optional{MaskGrid{NDims}}
    rule_masks::NTuple{NRules, Float32}

    applications::Vector{OrderedSet{CachedRuleApplication{NDims}}}
end

function RewriteCache(grid::CellGrid{NDims}, mask_grid::Optional{MaskGrid{NDims}},
                      rules::NTuple{NRules, Union{RewriteRule_Strip, RewriteRule_MD}},
                      mask_for_each_rule::NTuple{NRules, Float32},
                      context::MarkovOpContext
                     ) where {NDims, NRules}
    whole_grid_range = Box(min=one(CellIdx{NDims}), max=vsize(grid))

    @logic_logln("Caching valid rule applications...")
    @logic_tab_in()
    rule_info_tuple = map(rules, mask_for_each_rule) do rule, mask
        if rule isa RewriteRule_MD
            @bp_check(ndims(rule) <= NDims,
                      "Using ", ndims(rule), "D rule within a mere ", NDims, "D grid")
        end
        set = markov_allocator_acquire_ordered_set(context.allocator, CachedRuleApplication{NDims})
        empty!(set)

        @logic_logln("Masked ", mask, ", rule ", rules)
        @logic_tab_in()
        function process_rule_match(cell::CellIdx{NDims}, dir, is_matching::Bool)
            if is_matching
                @logic_log("   ", cell, " along ",
                           if dir isa GridDir
                               (dir.axis, "|", (dir.sign > 0 ? "+" : ""), dir.sign)
                           else
                               ("<", iter_join(dir.grid_to_rule, ", ")..., ">")
                           end...)
                @md_symm_logln() # If also logging MD-stuff then we need a newline between our match logs
                push!(set, (cell, dir))
            end
            return nothing
        end
        # Gather all matches whether this rule is a Strip or MD,
        #    and as a side effect grab the MD cached-orientations.
        md_orientations = if rule isa RewriteRule_Strip
            visit_rule_match_data(process_rule_match, rule,
                                  grid, mask_grid, mask, whole_grid_range)
            nothing
        else
            cached_orientations = RewriteRule_MD_Orientations(
                rule.cells, rule.symmetry,
                Val(NDims), context.allocator
            )
            visit_rule_match_data(process_rule_match, rule, cached_orientations,
                                  grid, mask_grid, mask, whole_grid_range)
            cached_orientations
        end
        @logic_logln()
        @logic_tab_out()

        return set, md_orientations
    end
    applications_tuple = Tuple(t[1] for t in rule_info_tuple)
    md_orientations_tuple = Tuple(t[2] for t in rule_info_tuple)
    @logic_logln("Total candidates: ", sum(map(length, applications_tuple), init=0))
    @logic_logln("MD orientations per rule: ", map(o -> (isnothing(o) ? "N/A" : o.n_entries), md_orientations_tuple))
    @logic_tab_out()

    applications = markov_allocator_acquire_array(
        context.allocator,
        tuple(length(applications_tuple)),
        OrderedSet{CachedRuleApplication{NDims}}
    )
    empty!(applications)
    append!(applications, applications_tuple)

    return RewriteCache{NDims, NRules, typeof(grid), typeof(rules)}(
        grid, rules, md_orientations_tuple,
        mask_grid, mask_for_each_rule,
        applications
    )
end

"Updates a rewrite cache, given an area of the grid that potentially changed"
function update_rewrite_cache!(cache::RewriteCache{NDims, NRules, TGrid, TRules},
                               range_to_invalidate::BoxI{NDims}
                              )::Nothing where {NDims, NRules, TGrid, TRules}
    foreach(cache.rules, cache.rule_masks, cache.applications, cache.rule_md_orientations
           ) do rule, mask, set, md_orientation_set
    #begin
        function process_match(cell, dir, is_matching::Bool)
            used_to_match = ((cell, dir) in set) #TODO: Improve equality/hashing performance for Orientation by adding an index counter
            if is_matching && !used_to_match
                push!(set, (cell, dir))
            elseif !is_matching && used_to_match
                delete!(set, (cell, dir))
            end
            return nothing
        end
        if rule isa RewriteRule_Strip
            visit_rule_match_data(process_match, rule,
                                  cache.grid, cache.mask_grid, mask, range_to_invalidate)
        else
            visit_rule_match_data(process_match,
                                  rule::RewriteRule_MD{NDims},
                                  md_orientation_set::RewriteRule_MD_Orientations{NDims},
                                  cache.grid, cache.mask_grid, mask, range_to_invalidate)
        end
        return nothing
    end
    return nothing
end

"Releases this cache's allocations to the allocator; don't touch this cache afterwards"
function close_rewrite_cache(cache::RewriteCache, context::MarkovOpContext)
    foreach(app -> markov_allocator_release_ordered_set(context.allocator, app), cache.applications)
    markov_allocator_release_array(context.allocator, cache.applications)
    foreach(o -> exists(o) && close(o), cache.rule_md_orientations)
    return nothing
end


#####################################
#  Rewrite Op: priorities interface

"A potential application of a specific rule to a specific part of the grid"
struct RewritePotentialApplication{NDims}
    rule_idx::Int
    desirability::Float32 # Weight*biases. May have any non-negative value.
    start_cell::CellIdx{NDims}
    dir::Union{GridDir, RewriteRule_MD_Orientation{NDims}}
end

"
Collective information about a group of potential rule applications.
Desirability is a potential application's bias value multiplied by the rule's weight.
"
struct RewriteGroupDesirability
    min::Float32
    max::Float32
    sum::Float32
end
RewriteGroupDesirability() = RewriteGroupDesirability(Inf32, -Inf32, 0.0f0)
function add_new_desirability(data::RewriteGroupDesirability, new_val::Float32)
    @markovjunior_assert(new_val >= 0, "Negative desirability: ", new_val, " among ", data)
    RewriteGroupDesirability(
        min(data.min, new_val),
        max(data.max, new_val),
        data.sum + new_val
    )
end
#

#TODO: Priorities should have control over evaluation of desirability, as most of them don't need to know the desirability of EVERY application for EVERY rule. This means they need to have an associated state object.

"
Decides which rewrite rule to apply.

Must implement the following interface:
* `pick_rule_using_rewrite_priority(...)` to return the rule index.
* `parse_markovjunior_rewrite_priority(::Val{:x}, ...)`
  to parse itself from the statement `PRIORITIZE(x, args...)`
* `dsl_string(self)` to turn the struct back into a DSL statement.
"
abstract type AbstractMarkovRewritePriority end
Base.:(==)(a::AbstractMarkovRewritePriority, b::AbstractMarkovRewritePriority) = (
    typeof(a) == typeof(b) &&
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))
)

pick_rule_using_rewrite_priority(priority::AbstractMarkovRewritePriority,
                                 rewite_op, rewrite_op_state,
                                 rng, op_context
                                )::Int = error(
    "Unimplemented: ", typeof(priority)
)
parse_markovjunior_rewrite_priority(::Val{Name}, expr_args, inputs::MacroParserInputs) where {Name} =
    error("Unimplemented: ", Name)

dsl_string(p::AbstractMarkovRewritePriority) = error("Unimplemented: ", typeof(p))


##################
#  Rewrite Op

"The beating heart of MarkovJunior: matches groups of pixels against a pattern, and rewrites them"
struct MarkovOpRewrite{TRules <: Tuple{Vararg{Union{RewriteRule_Strip, RewriteRule_MD}}},
                       TBias <: Tuple{Vararg{AbstractMarkovBias}},
                       TPriority <: AbstractMarkovRewritePriority
                      } <: AbstractMarkovOp
    priority::TPriority
    rules::TRules
    threshold::Optional{Threshold}
    biases::TBias
end

mutable struct MarkovOpRewrite_State{NDims, NRules,
                                     TGrid, TRules<:NTuple{NRules, Union{RewriteRule_Strip, RewriteRule_MD}},
                                     NBiases, TFullBias<:NTuple{NBiases, AbstractMarkovBias},
                                              TBiasStates<:NTuple{NBiases, Any}}
    # If not given, it can apply infinitely-many times.
    applications_left::Optional{Int}

    rewrite_cache::RewriteCache{NDims, NRules, TGrid, TRules}

    # All biases, incluing those inherited from parent ops.
    biases::TFullBias
    bias_states::TBiasStates

    #TODO: Rename the below fields to be more correct.

    # Every possible rule application: (rule_idx, priority, first_cell, dir).
    # They are sorted by rule -- see 'weighted_options_buffer_first_indices'.
    weighted_options_buffer::Vector{RewritePotentialApplication{NDims}}
    # For each rule, indicates the first index of the first entry
    #    in 'weighted_options_buffer' using that rule.
    #
    # An extra entry is at the end to simplify use of this array.
    # The range of entries for rule i (abbreviating fields as 'wobfi' and 'wob')
    #    is 'wobfi[i] : (wobfi[i+1]-1)'.
    weighted_options_buffer_first_indices::Vector{Int}

    # Info about the range of rule weights*biases, across all rules.
    weight_data_buffer::RewriteGroupDesirability
    # Info about all potential applications of each individual rule.
    weight_data_buffer_per_rule::Vector{RewriteGroupDesirability}

    # Buffer for holding a final group of sub-choices to choose from.
    final_choices_buffer::Vector{RewritePotentialApplication{NDims}}
    # Buffer for holding the previous pixel values after applying a 1D rewrite rule.
    # At least as long as the longest 1D rewrite rule.
    prev_1D_pixels_buffer::Vector{UInt8}
end
rewrite_rule_option_indices(state::MarkovOpRewrite_State, rule_idx::Integer) = (
    state.weighted_options_buffer_first_indices[rule_idx] :
     (state.weighted_options_buffer_first_indices[rule_idx+1] - 1)
)

function markov_op_state_type(op::MarkovOpRewrite{TRules}, TGrid::Type{<:CellGrid{NGridDims}},
                              rng::PRNG, context::MarkovOpContext
                             ) where {NGridDims, TRules}
    all_biases = Iterators.flatten((context.all_biases, op.biases))
    return MarkovOpRewrite_State{
        NGridDims, length(op.rules), TGrid, TRules,
        count(i->true, all_biases), Tuple{typeof.(all_biases)...},
        Tuple{markov_bias_state_type.(all_biases, Ref(TGrid), Ref(rng), Ref(context.bias_context))...}
    }
end
function markov_op_initialize(r::MarkovOpRewrite{<:NTuple{NRules, Any}, TBias, TPriority},
                              TState::Type{<:MarkovOpRewrite_State{NDims, NRules, TGrid, TRules,
                                                                   NFullBiases, TFullBiases, TFullBiasStates}},
                              grid::TGrid, rng::PRNG,
                              context::MarkovOpContext
                             ) where {NDims, NRules, TGrid<:CellGrid{NDims}, TRules, TBias, TPriority,
                                      NFullBiases, TFullBiases, TFullBiasStates}
    mask_grid = if all(rule -> isnothing(rule.mask), r.rules)
        nothing
    else
        a = markov_allocator_acquire_array(context.allocator, size(grid), Float32)
        rand!(rng, a)
        a
    end

    cache = RewriteCache(grid, mask_grid, r.rules,
                         map(ru -> pick_mask(ru, rng), r.rules),
                         context)

    append!(context.all_biases, r.biases)
    biases::TFullBiases = Tuple(context.all_biases)
    bias_states::TFullBiasStates = Tuple(Iterators.map(
        (b,t) -> markov_bias_initialize(b, t, grid, rng, context.bias_context),
        biases, TFullBiasStates.parameters
    ))

    threshold = isnothing(r.threshold) ? nothing : get_threshold(r.threshold, grid, rng)
    longest_1D_rule::Int = maximum(Tuple(
        length(ru.cells) for ru in r.rules if ru isa RewriteRule_Strip
    ); init=0)

    out_state = TState(
        threshold, cache,
        biases, bias_states,
        markov_allocator_acquire_array(context.allocator, tuple(128), RewritePotentialApplication{NDims}),
        markov_allocator_acquire_array(context.allocator, tuple(NRules+1), Int),
        RewriteGroupDesirability(),
        markov_allocator_acquire_array(context.allocator, tuple(NRules), RewriteGroupDesirability),
        markov_allocator_acquire_array(context.allocator, tuple(256), RewritePotentialApplication{NDims}),
        markov_allocator_acquire_array(context.allocator, tuple(max(longest_1D_rule, 1)), UInt8)
    )
    if all(isempty, cache.applications)
        @logic_logln("MarkovOpRewrite has no options at the start; canceling...")
        markov_op_cancel(r, out_state, context)
        return nothing
    else
        @logic_logln("MarkovOpRewrite has been initialized; there are ",
                      sum(map(length, cache.applications), init=0), " initial options")
        @logic_logln("The actual threshold is ",
                      typeof(out_state.applications_left), "(", out_state.applications_left, ")")
        return out_state
    end
end
function markov_op_iterate(r::MarkovOpRewrite{TRules, TSelfBiases, TPriority},
                           state::MarkovOpRewrite_State{NDims, NRules, TGrid, TRules, NBiases, TFullBias, TBiasStates},
                           grid::CellGrid{NDims}, rng::PRNG,
                           context::MarkovOpContext,
                           ticks_left::Ref{Optional{Int}}
                          ) where {NDims, NRules, TGrid, TRules, TPriority,
                                   NBiases, TFullBias, TBiasStates, TSelfBiases}
    @markovjunior_assert(get_something(state.applications_left, 1) > 0)
    @logic_logln("Remaining threshold for this rewrite op: ",
                get_something(state.applications_left, "infinity"))

    while isnothing(ticks_left[]) || (ticks_left[] > 0)
        # Get all the options and their weights.
        empty!(state.weighted_options_buffer)
        state.weighted_options_buffer_first_indices .= -1
        state.weight_data_buffer = RewriteGroupDesirability()
        for i in 1:NRules
            state.weight_data_buffer_per_rule[i] = RewriteGroupDesirability()
        end
        foreach(r.rules, 1:NRules) do rule, rule_i
            state.weighted_options_buffer_first_indices[rule_i] = 1 + length(state.weighted_options_buffer)
            for (start_cell, dir) in state.rewrite_cache.applications[rule_i]
                cell_at = if rule isa RewriteRule_Strip
                    CellLine(start_cell, dir, convert(Int32, length(rule.cells)))
                elseif rule isa RewriteRule_MD
                    rule_size = vsize((dir::RewriteRule_MD_Orientation{NDims}).rule_permutation)
                    CellRegion(start_cell:(start_cell + rule_size))
                else
                    error("Unhandled: ", typeof(rule))
                end
                biases::NTuple{NBiases, Optional{Float32}} = markov_bias_calculate.(
                    state.biases, state.bias_states,
                    Ref(grid), Ref(cell_at), Ref(rng)
                )
                @markovjunior_assert(all(b -> something(b, 1.0f0) >= 0, biases),
                                     "Some biases returned negative values! ",
                                       collect(zip(typeof.(state.biases), biases)))

                if all(exists, biases)
                    desirability = rule.weight * (iszero(NBiases) ? 1.0f0 : sum(biases, init=zero(Float32)))
                    push!(state.weighted_options_buffer, RewritePotentialApplication(
                        rule_i, desirability,
                        start_cell, dir
                    ))

                    state.weight_data_buffer = add_new_desirability(
                        state.weight_data_buffer, desirability
                    )
                    state.weight_data_buffer_per_rule[rule_i] = add_new_desirability(
                        state.weight_data_buffer_per_rule[rule_i], desirability
                    )
                end
            end
        end
        state.weighted_options_buffer_first_indices[NRules+1] = 1 + length(state.weighted_options_buffer)

        @logic_logln("There are ", length(state.weighted_options_buffer),
                      " options with biases ranging from ",
                      state.weight_data_buffer.min, " to ", state.weight_data_buffer.max,
                      ": [")
        log_logic() && for option::RewritePotentialApplication{NDims} in state.weighted_options_buffer
            @logic_logln("\tRule ", option.rule_idx, " at ", option.start_cell, " along ",
                         (option.dir isa GridDir) ?
                             option.dir :
                             Tuple((g.axis * g.sign) for g in option.dir.rule_to_grid),
                         "; desirability=", option.desirability)
        end
        @logic_logln("]\n")

        # If no options are left, we're done.
        if isempty(state.weighted_options_buffer)
            markov_op_cancel(r, state, context)
            return nothing
        end

        # Pick a rule using the chosen Priority.
        # Interpret an invalid choice as "no options left".
        pick_rule_i = pick_rule_using_rewrite_priority(r.priority, r, state, rng, context)
        @logic_logln("Chose rule ", pick_rule_i)
        if !in(pick_rule_i, 1:NRules)
            @logic_logln("Invalid rule index! This is the end of the Op.")
            markov_op_cancel(r, state, context)
            return nothing
        end
        pick_options_range = rewrite_rule_option_indices(state, pick_rule_i)
        if isempty(pick_options_range)
            @logic_logln("Rule has no options! This is the end of the Op.")
            markov_op_cancel(r, state, context)
            return nothing
        end
        picked_options = @view state.weighted_options_buffer[pick_options_range]
        @markovjunior_assert(all(o -> o.rule_idx == pick_rule_i, picked_options),
                             "We chose rule ", pick_rule_i, " but its options (",
                               pick_options_range, ") include rules ",
                               unique(o->o.rule_idx for o in picked_options))
        rule_desirability_data = state.weight_data_buffer_per_rule[pick_rule_i]

        # Pick an option within that rule.
        if rule_desirability_data.min == rule_desirability_data.max
            @logic_logln("No biases exist, so we'll use uniform-random selection")
        else
            @logic_logln("Picking the move with the highest desirability (",
                         rule_desirability_data.max, ")")
            empty!(state.final_choices_buffer)
            append!(state.final_choices_buffer,
                    Iterators.filter(app -> app.desirability == rule_desirability_data.max,
                                     picked_options))
            picked_options = @view state.final_choices_buffer[1:end] # Explicit indices for type-stability
            @markovjunior_assert(!isempty(picked_options),
                                 "Somehow none of the options in this rule have their max desirability ",
                                   rule_desirability_data.max)
            if length(picked_options) > 1
                @logic_logln("\tThere are ", length(picked_options), ", so one will be chosen randomly")
            end
        end
        (pick_start_cell, pick_dir) = let p = rand(rng, picked_options)
            (p.start_cell, p.dir)
        end
        @logic_logln("Decided to apply the rule at ", pick_start_cell, " along ",
                     if pick_dir isa GridDir
                        ((pick_dir.sign < 0 ? '-' : '+'), pick_dir.axis)
                     else
                        ("<", iter_join(pick_dir.rule_to_grid, ", ")..., ">")
                     end...)

        # Apply the rule.
        # Because each rule is a different type but known at compile-time,
        #    we should add a layer of dispatch when executing it.
        (rule -> begin
            (affected_area, original_values) = if rule isa RewriteRule_Strip
                source_values = Tuple(
                    grid[grid_dir_pos_along(pick_dir, pick_start_cell, i-1)]
                    for i in 1:length(rule.cells)
                )
                # Each rule's rewrite cell is also a different type known at compile-time.
                foreach(rule.cells, 1:length(rule.cells)) do (rewrite_source, rewrite_dest), cell_i
                    cell_pos = grid_dir_pos_along(pick_dir, pick_start_cell, cell_i-1)
                    state.prev_1D_pixels_buffer[cell_i] = grid[cell_pos]
                    grid[cell_pos] = pick_rewrite_value(rewrite_dest, rewrite_source,
                                                        source_values, cell_i,
                                                        rng)
                end

                # Calculate the affected area.
                rule_len = convert(Int32, length(rule.cells))
                pick_end_cell = grid_dir_pos_along(pick_dir, pick_start_cell, rule_len-1)
                (pick_start_cell, pick_end_cell) = minmax(pick_start_cell, pick_end_cell)
                b = Box(pick_start_cell:pick_end_cell)

                # Create an N-dimensional view of the source values in the grid.
                source_view_shape = ntuple(Val(NDims)) do i
                    if i == pick_dir.axis
                        Base.:(:)
                    else
                        1
                    end
                end
                source_view = reshape(@view(state.prev_1D_pixels_buffer[1:end]),
                                      source_view_shape...)

                b, source_view
            else
                pick_dir_o::RewriteRule_MD_Orientation{NDims} = pick_dir

                # Grab a view of the grid area covered by the rule.
                values_range = ntuple(Val(NDims)) do i
                    pick_start_cell[i] : (pick_start_cell[i] + size(pick_dir_o.rule_permutation, i) - 1)
                end
                values_view = @view grid[values_range...]

                # Grab a buffer of the grid area's original values.
                source_values = pick_dir_o.application_buffer
                source_values .= values_view

                # Generate the grid area's new values.
                values_view .= pick_rewrite_value.((t->t[2]).(pick_dir_o.rule_permutation),
                                                   (t->t[1]).(pick_dir_o.rule_permutation),
                                                   Ref(source_values),
                                                   Tuple.(eachindex(IndexCartesian(), pick_dir_o.rule_permutation)),
                                                   Ref(rng))

                # Calculate the affected area.
                BoxI{NDims}(
                    min=pick_start_cell,
                    size=vsize(pick_dir.rule_permutation)
                ), source_values
            end

            # Update inner bookkeeping.
            update_rewrite_cache!(state.rewrite_cache, affected_area)
            state.bias_states = markov_bias_update.(
                state.biases, state.bias_states,
                Ref(grid), Ref(affected_area), Ref(original_values),
                Ref(rng)
            )
        end)(r.rules[pick_rule_i])

        # Update outer bookkeeping.
        if exists(ticks_left[])
            ticks_left[] -= 1
        end
        if exists(state.applications_left)
            state.applications_left -= 1
            if state.applications_left < 1
                markov_op_cancel(r, state, context)
                return nothing
            end
        end
    end

    return state
end

function markov_op_cancel(op::MarkovOpRewrite, s::MarkovOpRewrite_State,
                          context::MarkovOpContext)
    # We're freeing a lot of stuff here, so move the allocator to a type-stable context.
    function run(allocator::T) where {T<:AbstractMarkovAllocator}
        foreach(markov_bias_cleanup, s.biases, s.bias_states)
        close_rewrite_cache(s.rewrite_cache, context)

        if exists(s.rewrite_cache.mask_grid)
            markov_allocator_release_array(allocator, s.rewrite_cache.mask_grid)
        end
        markov_allocator_release_array(allocator, s.weighted_options_buffer)
        markov_allocator_release_array(allocator, s.weight_data_buffer_per_rule)
        markov_allocator_release_array(allocator, s.weighted_options_buffer_first_indices)
        markov_allocator_release_array(allocator, s.final_choices_buffer)
        markov_allocator_release_array(allocator, s.prev_1D_pixels_buffer)
    end
    run(context.allocator)
end


#####################
#  DSL integration

dsl_string_rewrite_source(rule::UInt8) = dsl_string(rule)
dsl_string_rewrite_source(rule::RewriteRuleCell_Set) = "[$(dsl_string(rule))]"
dsl_string_rewrite_source(rule::RewriteRuleCell_Wildcard) = "_"

dsl_string_rewrite_dest(rule::UInt8) = dsl_string(rule)
dsl_string_rewrite_dest(rule::RewriteRuleCell_Set) = "{$(dsl_string(rule))}"
dsl_string_rewrite_dest(rule::RewriteRuleCell_List) = "[$(string(dsl_string.(rule)...))]"
dsl_string_rewrite_dest(rule::RewriteRuleCell_Wildcard) = "_"
dsl_string_rewrite_dest(rule::RewriteRuleCell_Lookup{Int}) = "[$(rule.source_idx)]"
dsl_string_rewrite_dest(rule::RewriteRuleCell_Lookup{<:NTuple}) = "[$(iter_join(rule.source_idx, ", ")...)]"

dsl_string_rewrite_mask(mask::Nothing) = ""
dsl_string_rewrite_mask(mask::Float32) = "%$mask"
dsl_string_rewrite_mask(mask::NTuple{2, Float32}) = "%($(mask[1]):$(mask[2]))"

function dsl_string_rewrite_md_array(array::Array{RewriteCell_MD{NDims}, NDims},
                                     take_source::Bool) where {NDims}
    output = preallocated_vector(Char, 1024)
    append!(output, "[\n  ")

    # Define the incrementing index counter to iterate through the array in order
    #    and track which axis boundary we just crossed.
    idx = MVector(ntuple(one, Val(NDims)))
    function advance_idx_and_get_changed_axis()::Int
        for i in 1:NDims
            idx[i] += 1
            if idx[i] <= size(array, i)
                return i
            else
                idx[i] = 1
            end
        end

        return NDims+1
    end

    # Iterate and generate.
    changed_axis = 1
    while changed_axis <= NDims
        # Write the separator.
        # Note that the first axis from the user's POV is the column, not the row.
        if changed_axis == 1
            append!(output, " ")
        elseif changed_axis == 2
            append!(output, "\n  ")
        else
            push!(output, ' ')
            for _ in 1:changed_axis
                push!(output, ';')
            end
            append!(output, "\n  ")
        end

        # Write the current element.
        append!(output,
            if take_source
                dsl_string_rewrite_source(array[idx...][1])
            else
                dsl_string_rewrite_dest(array[idx...][2])
            end
        )

        changed_axis = advance_idx_and_get_changed_axis()
    end

    append!(output, "\n ]")
    return String(output)
end

dsl_string_rewrite_axis(a::Integer) = if a < 5
    ('x', 'y', 'z', 'w')[a]
else
    a
end
function dsl_string_rewrite_signed_axis(a::Integer)
    letter = if abs(a) in 1:4
        ('x', 'y', 'z', 'w')[abs(a)]
    else
        nothing
    end
    return string(
        if a > 0
            '+'
        else
            '-'
        end,
        if isnothing(letter)
            "("
        else
            ""
        end,
        if exists(letter)
            letter
        else
            abs(a)
        end,
        if isnothing(letter)
            ")"
        else
            ""
        end
    )
end
dsl_string(@nospecialize strip::RewriteRule_Strip) = string(
    dsl_string_rewrite_source.(t[1] for t in strip.cells)...,
    " => ",
    dsl_string_rewrite_dest.(t[2] for t in strip.cells)...,
    " $(dsl_string_rewrite_mask(strip.mask))",
    (isone(strip.weight) ? () : (" *", strip.weight))...,
    if isempty(strip.explicit_symmetries) && isnothing(strip.tail_symmetry)
        ()
    else
        (
            " \\[ ",
            # Explicit symmetries:
            iter_join_flatten(
                (dsl_string_rewrite_signed_axis(dir.axis * dir.sign)
                  for dir in strip.explicit_symmetries),
                ", "
            )...,
            # Unlimited symmetries:
            if exists(strip.tail_symmetry) && !isempty(strip.explicit_symmetries)
                ", "
            else
                ""
            end,
            if isnothing(strip.tail_symmetry)
                tuple()
            elseif strip.tail_symmetry isa Int
                tuple(
                    strip.tail_symmetry, "..."
                )
            elseif strip.tail_symmetry isa NTuple{2, Optional{Int}}
                tuple(
                    if exists(strip.tail_symmetry[1])
                        ("-(", strip.tail_symmetry[1], ")...",
                         exists(strip.tail_symmetry[2]) ? ", " : "")
                    else
                        ()
                    end...,
                    if exists(strip.tail_symmetry[2])
                        ("+(", strip.tail_symmetry[2], ")...")
                    else
                        ()
                    end...
                )
            else
                error("Unhandled: ", typeof(strip.tail_symmetry))
            end...,
            " ]"
        )
    end...
)
dsl_string(@nospecialize md::RewriteRule_MD) = string(
    dsl_string_rewrite_md_array(md.cells, true),
    " => ",
    dsl_string_rewrite_md_array(md.cells, false),
    " $(dsl_string_rewrite_mask(md.mask))",
    (isone(md.weight) ? () : (" *", md.weight))...,
    if isempty(md.symmetry.grid_axis_choices) && isempty(md.symmetry.chiral_groups)
        ()
    else
        (
            " \\[\n",

            # Chiral groups:
            iter_join_flatten(
                Iterators.map(md.symmetry.chiral_groups) do cg
                    return tuple(
                        '{',
                        iter_join(
                            Iterators.map(dsl_string_rewrite_axis, cg),
                            ", "
                        )...,
                        '}'
                    )
                end,
                ",\n    "
            )...,

            # Delimeter between chiral groups and axis choices:
            if !isempty(md.symmetry.grid_axis_choices) && !isempty(md.symmetry.chiral_groups)
                ",\n    "
            else
                ""
            end,

            # Axis choices:
            iter_join_flatten(
                Iterators.map(md.symmetry.grid_axis_choices) do (perm_data, tail_symmetry)
                    rule_axes = @view(perm_data[:, 1])
                    permutations = @view(perm_data[:, 2:end])
                    is_multi = (length(rule_axes) > 1)
                    return tuple(
                        # Header:
                        if is_multi
                            tuple(
                                '(',
                                iter_join(
                                    Iterators.map(dsl_string_rewrite_axis, rule_axes),
                                    ", "
                                )...,
                                ')'
                            )
                        else
                            tuple(
                                dsl_string_rewrite_axis(rule_axes[1])
                            )
                        end...,
                        "[ ",

                        # Elements:
                        if is_multi
                            iter_join_flatten(
                                Iterators.map(1:size(permutations, 2)) do perm_i
                                    permutation = @view permutations[:, perm_i]
                                    return tuple(
                                        '(',
                                        iter_join(
                                            Iterators.map(permutation) do grid_dir_int
                                                return dsl_string_rewrite_signed_axis(grid_dir_int)
                                            end,
                                            ", "
                                        )...,
                                        ')'
                                    )
                                end,
                                ", "
                            )
                        else
                            iter_join(
                                Iterators.map(dsl_string_rewrite_signed_axis, permutations),
                                ", "
                            )
                        end...,

                        # Separator between elements and tail symmetry:
                        if !isempty(permutations) && !isa(tail_symmetry, Union{Nothing, NTuple{2, Nothing}})
                            ", "
                        else
                            ""
                        end,

                        # Tail symmetries:
                        if isnothing(tail_symmetry)
                            tuple()
                        elseif tail_symmetry isa Int
                            tuple(
                                tail_symmetry, "..."
                            )
                        elseif tail_symmetry isa NTuple{2, Optional{Int}}
                            tuple(
                                if exists(tail_symmetry[1])
                                    ("-(", tail_symmetry[1], ")...",
                                     exists(tail_symmetry[2]) ? ", " : "")
                                else
                                    ()
                                end...,
                                if exists(tail_symmetry[2])
                                    ("+(", tail_symmetry[2], ")...")
                                else
                                    ()
                                end...
                            )
                        else
                            error("Unhandled: ", typeof(tail_symmetry))
                        end...,

                        " ]"
                    )
                end,
                ",\n    "
            )...,

            "\n]"
        )
    end...
)

dsl_string(@nospecialize op::MarkovOpRewrite) = string(
    "@rewrite ",
    exists(op.threshold) ? dsl_string(op.threshold) : "",
      " ",
    if length(op.rules) == 1
        dsl_string(op.rules[1])
    else
        string(
            "begin\n    ",
            "PRIORITIZE(", dsl_string(op.priority), ")\n    ",
            iter_join(dsl_string.(op.rules), "\n    ")...,
            "\nend"
        )
    end,
      " ",
    if length(op.biases) == 1
        dsl_string(op.biases[1])
    elseif length(op.biases) > 1
        string(
            "begin\n    ",
            iter_join(dsl_string.(op.biases), "\n    ")...,
            "\nend"
        )
    else
        ""
    end
)

"Turns an axis name/index into its index"
function parse_markovjunior_axis_name(loc, inputs::MacroParserInputs, expr)::Int
    if expr == :x
        return 1
    elseif expr == :y
        return 2
    elseif expr == :z
        return 3
    elseif expr == :w
        return 4
    elseif !isa(expr, Integer)
        raise_parse_error(loc, inputs, "Axis must be a number of one of [ x, y, z, w ]! Got `", expr, "`")
    elseif expr < 1
        raise_parse_error(loc, inputs, "Axis must be a 1-based index! Got ", expr)
    else
        return convert(Int, expr)
    end
end
"
Turns a signed axis name (`x`, `-(1)`, `4`, `+z`, etc) into a tuple of axis and sign.
If no sign is provided, the sign is set to 0.

You can optionally have this return `nothing` when invalid syntax is given (by default it throws an error).
"
function parse_markovjunior_axis_signed(loc, inputs::MacroParserInputs, expr,
                                        error_as_fallback::Bool = true)::Optional{NTuple{2, Int}}
    if (expr in (:x, :y, :z, :w)) || (expr isa Integer && expr > 0)
        return (parse_markovjunior_axis_name(loc, inputs, expr), 0)
    elseif (expr isa Integer) && iszero(expr)
        raise_parse_error(loc, inputs, "Axis must be a 1-based value!")
    elseif (expr isa Integer) # Must be negative
        raise_parse_error(loc, inputs, "Signed axis syntax requires parentheses, e.g. `+(3)` instead of `+3`")
    elseif @capture(expr, -(a_Integer))
        return (a, -1)
    elseif @capture(expr, +(a_Integer))
        return (a, +1)
    elseif @capture(expr, -(x))
        return (1, -1)
    elseif @capture(expr, -(y))
        return (2, -1)
    elseif @capture(expr, -(z))
        return (3, -1)
    elseif @capture(expr, -(w))
        return (4, -1)
    elseif @capture(expr, +(x))
        return (1, +1)
    elseif @capture(expr, +(y))
        return (2, +1)
    elseif @capture(expr, +(z))
        return (3, +1)
    elseif @capture(expr, +(w))
        return (4, +1)
    elseif error_as_fallback
        raise_parse_error(loc, inputs, "Unexpected syntax: `", expr, "`")
    else
        nothing
    end
end

"Note that Source (lhs) patterns will return a Set as a List, to preserve input order for later processing"
function parse_markovjunior_rewrite_rule_strip_side(inputs::MacroParserInputs, loc, expr,
                                                    is_source::Bool)::Vector
    push!(inputs.op_stack_trace, is_source ? "Left-hand side" : "Right-hand side")
    try
        try_lookup_char(c::Char)::Union{RewriteRuleCellSource, RewriteRuleCellDest{Int}} = if haskey(CELL_CODE_BY_CHAR, c)
            CELL_CODE_BY_CHAR[c]
        elseif c == '_'
            RewriteRuleCell_Wildcard()
        else
            raise_parse_error(
                loc, inputs,
                "Unsupported color value '", c, "'! Supported are '_' and ",
                  "[ ", iter_join(keys(CELL_CODE_BY_CHAR), ", ")..., " ]"
            )
        end

        # The complex sequence of pixel rules ultimately turns into a series of binary operators --
        #   implicit multiplication, array indexing, braces.
        # Start by flattening this AST into a simpler representation:
        #     a list of Char, :set=>Vector{Char}, :list=>Vector{Char}, and :ref=>Int.
        flatten_rule_expr(s::Symbol) = collect(string(s))
        function flatten_rule_expr(e::Expr)
            if @capture(e, a_{b_})
                [ flatten_rule_expr(a)..., :set=>flatten_rule_expr(b) ]
            elseif @capture(e, a_[b_])
                [
                    flatten_rule_expr(a)...,
                    if b isa Integer
                        :ref=>convert(Int, b)
                    else
                        :list=>flatten_rule_expr(b)
                    end
                ]
            elseif @capture(e, [ a_ ])
                if a isa Integer
                    [ :ref=>convert(Int, a) ]
                else
                    [ :list=>flatten_rule_expr(a) ]
                end
            elseif @capture(e, {a_})
                [ :set=>flatten_rule_expr(a) ]
            elseif @capture(e, a_*b_)
                [ flatten_rule_expr(a)..., flatten_rule_expr(b)... ]
            else
                raise_parse_error(loc, inputs, "Unsupported bit of syntax: ", e)
            end
        end
        flatten_rule_expr(o) = raise_parse_error(loc, inputs, "Unsupported expression: ", typeof(o), "(", o, ")")

        flattened = flatten_rule_expr(expr)

        # Now turn each element into a proper data representation.
        return collect(Union{RewriteRuleCellSource, RewriteRuleCellDest{Int}},
                       Iterators.map(flattened) do simple_repr
            if simple_repr isa Char
                return try_lookup_char(simple_repr)
            # If processing the Source side, we should keep sets as lists for now
            #   as their order may be relevant to the Destination side.
            elseif (simple_repr isa Pair) && (simple_repr[1] == :list)
                invalid_elements = filter(c -> !isa(c, Char), simple_repr[2])
                if isempty(invalid_elements)
                    return RewriteRuleCell_List(Tuple(try_lookup_char.(simple_repr[2])))
                else
                    raise_parse_error(loc, inputs,
                                   "Invalid nested syntax in ", is_source ? "Source" : "Dest",
                                     " part of rule!")
                end
            elseif (simple_repr isa Pair) && (simple_repr[1] == :set)
                if is_source
                    raise_parse_error(loc, inputs, "'Set' syntax (`{RGB}`) not allowed on the Source side")
                end

                invalid_elements = filter(c -> !isa(c, Char), simple_repr[2])
                if isempty(invalid_elements)
                    return RewriteRuleCell_Set(simple_repr[2]...)
                else
                    raise_parse_error(loc, inputs,
                                   "Invalid nested syntax in ", is_source ? "Source" : "Dest",
                                     " part of rule!")
                end
            elseif (simple_repr isa Pair) && (simple_repr[1] == :ref)
                if !isa(simple_repr[2], Integer)
                    raise_parse_error(loc, inputs,
                                      "Expected source-ref to look like `[5]`, got `[", simple_repr[2], "]`")
                else
                    return RewriteRuleCell_Lookup(convert(Int, simple_repr[2]))
                end
            else
                @bp_check(false, "Unexpected: ", simple_repr)
            end
        end)
    finally
        pop!(inputs.op_stack_trace)
    end
end
function parse_markovjunior_rewrite_rule_strip_symmetry(inputs::MacroParserInputs, loc, exprs)
    # If no symmetry is given at all, then support all symmetries.
    if isnothing(exprs)
        return (GridDir[ ], 1)
    end

    push!(inputs.op_stack_trace, "symmetry statement")
    try
        s_explicit = GridDir[ ]
        s_tail::RewriteRule_TailSymmetry = nothing
        get_axis(a) = parse_markovjunior_axis_name(loc, inputs, a)
        for expr in exprs
            if @capture(expr, a_...)
                # Get the axis and optional sign.
                (axis, dir) = parse_markovjunior_axis_signed(loc, inputs, a)

                # Incorporate the symmetry tail.
                if isnothing(s_tail)
                    if dir == 0
                        s_tail = axis
                    elseif dir == -1
                        s_tail = (axis, nothing)
                    else
                        @bp_check(dir == 1)
                        s_tail = (nothing, axis)
                    end
                elseif s_tail isa Int
                    raise_parse_error(loc, inputs,
                                      "Found both `", s_tail, "...` and `", expr,
                                        "` in one symmetry statement! They overlap each other")
                elseif (dir == -1)
                    if exists(s_tail[1])
                        raise_parse_error(loc, inputs,
                                          "Found both `-(", s_tail[1], ")...` and `",
                                            a, "` in one symmetry statement!")
                    else
                        s_tail = (axis, s_tail[2])
                    end
                elseif (dir == 1)
                    if exists(s_tail[2])
                        raise_parse_error(loc, inputs,
                                          "Found both `+(", s_tail[2], ")...` and `",
                                            a, "` in one symmetry statement!")
                    else
                        s_tail = (s_tail[1], axis)
                    end
                else
                    error("Unhandled: ", axis, "/", dir, "/", sprint(io -> dump(io, a, maxdepth=100)))
                end
            else
                (axis, dir) = parse_markovjunior_axis_signed(loc, inputs, expr)
                if dir in -1:0
                    push!(s_explicit, GridDir(axis, -1))
                end
                if dir in 0:1
                    push!(s_explicit, GridDir(axis, 1))
                end
            end
        end

        # Validate the parsed data.
        #   * No duplicates should exist:
        duplicates = Set{GridDir}()
        for (i, g) in enumerate(s_explicit)
            for j in (i+1):length(s_explicit)
                if g == s_explicit[j]
                    push!(duplicates, g)
                end
            end
        end
        if !isempty(duplicates)
            raise_parse_error(loc, inputs, "Duplicate symmetry values found: ", collect(duplicates))
        end
        #   * The lower-bound should start after the explicit list:
        tail_signs = ('-', '+')
        max_explicit_axis = (
            maximum((g.axis for g in s_explicit if g.sign<0), init=-1),
            maximum((g.axis for g in s_explicit if g.sign>0), init=-1)
        )
        tail_axes = if isnothing(s_tail)
            (nothing, nothing)
        elseif s_tail isa Int
            # The symmetry tail does not mention sign, so factor the sign out of the error-check.
            max_explicit_axis = max(max_explicit_axis...)
            tail_signs = tuple("")
            s_tail
        elseif s_tail isa NTuple{2, Optional{Int}}
            s_tail
        else
            error("Unhandled: ", typeof(s_tail))
        end
        if any(t -> exists(t[2]) && (t[1] >= t[2]), zip(max_explicit_axis, tail_axes))
            raise_parse_error(loc, inputs,
                              "Intersection between explicit symmetry axes (which go up to ",
                                max_explicit_axis, ") and tail symmetry axes (",
                                iter_join((string("`", sign, t, "...`") for (sign, t) in zip(tail_signs, tail_axes) if exists(t)),
                                          " and "))
        end

        return s_explicit, s_tail
    finally
        pop!(inputs.op_stack_trace)
    end
end
function parse_markovjunior_rewrite_rule_strip(inputs::MacroParserInputs, loc, expr)
    push!(inputs.op_stack_trace, "Rewrite rule `$expr`")
    try
        if !@capture(expr, lhsExpr_ => b_)
            raise_parse_error(loc, inputs, "Invalid format; expected `src => dest [modifiers]`")
        end

        # Strip out the modifiers.
        # NOTE: MacroTools has a bug where the | operator doesn't always work.
        weightMulScalar = nothing
        weightDivScalar = nothing
        symmetryExprs = nothing
        maskExpr = nothing
        @capture(    b, c_ % maskExpr_ * weightMulScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ / weightDivScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ \ [ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ * weightMulScalar_Real) ||
            @capture(b, c_ % maskExpr_ / weightDivScalar_Real) ||
            @capture(b, c_ * weightMulScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ / weightDivScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ \[ symmetryExprs__ ]) ||
            @capture(b, c_ * weightMulScalar_Real) ||
            @capture(b, c_ / weightDivScalar_Real) ||
            @capture(b, c_ % maskExpr_) ||
            (c = b)
        rhsExpr = c

        # Parse the modifiers.
        weight = convert(Float32, if exists(weightMulScalar)
            @bp_check(isnothing(weightDivScalar))
            weightMulScalar
        elseif exists(weightDivScalar)
            @bp_check(isnothing(weightMulScalar))
            1 / convert(Float64, weightDivScalar)
        else
            1
        end)
        mask = if @capture(maskExpr, maskAExpr_Real:maskBExpr_Real)
            convert.(Ref(Float32), (maskAExpr, maskBExpr))
        elseif maskExpr isa Real
            convert(Float32, maskExpr)
        elseif isnothing(maskExpr)
            nothing
        else
            raise_parse_error(loc, inputs, "Expected mask to be `%x` or `%(x:y)`; got `%$maskExpr`")
        end
        (symmetries_explicit, symmetries_tail) = parse_markovjunior_rewrite_rule_strip_symmetry(inputs, loc, symmetryExprs)

        # Parse the rule.
        lhs = parse_markovjunior_rewrite_rule_strip_side(inputs, loc, lhsExpr, true)
        rhs = parse_markovjunior_rewrite_rule_strip_side(inputs, loc, rhsExpr, false)

        # Post-process and validate the rule.
        if length(lhs) != length(rhs)
            raise_parse_error(loc, inputs,
                           "Source has ", length(lhs), " entries while Dest has ", length(rhs))
        end
        for i in 1:length(lhs)
            # If source cell is currently a List, remake it as a Set.
            # If dest cell is also a List, then its order needs to update accordingly.
            if lhs[i] isa RewriteRuleCell_List
                if rhs[i] isa RewriteRuleCell_List
                    if length(lhs[i]) != length(rhs[i])
                        raise_parse_error(loc, inputs,
                                       "Destination Cell ", i, " is a list of ", length(rhs[i]),
                                       " elements, but its Source cell has ", length(lhs[i]),
                                       " elements instead")
                    else
                        lhs_to_rhs = Dict{UInt8, UInt8}(
                            L => R for (L, R) in zip(lhs[i], rhs[i])
                        )
                        lhs[i] = RewriteRuleCell_Set(lhs[i]...)
                        rhs[i] = RewriteRuleCell_List(Tuple(lhs_to_rhs[L] for L in lhs[i]))
                    end
                else
                    lhs[i] = RewriteRuleCell_Set(lhs[i]...)
                end
            end

            # If dest cell is a list, the source cell must be too (at this point actually a Set).
            if (rhs[i] isa RewriteRuleCell_List) && !isa(lhs[i], RewriteRuleCell_Set)
                raise_parse_error(loc, inputs, "Destination Cell ", i, " is a list, but its Source cell is not!")
            end

            # If dest cell references a source cell, that reference must be valid.
            if (rhs[i] isa RewriteRuleCell_Lookup{Int}) && !in(rhs[i].source_idx, 1:length(lhs))
                raise_parse_error(loc, inputs,
                               "Destination Cell ", i, " references nonexistent source cell ",
                                 rhs[i].source_idx)
            end
        end

        # If the rule has only one cell, remove all symmetries to avoid redundancy.
        if length(rhs) == 1
            symmetries_explicit = [ GridDir(1, 1) ]
            symmetries_tail = nothing
        end
        #TODO: If the rule is symmetric (see requirements below), convert all symmetries to the positive side.
        #       * For source strip, is cell i identical to cell LEN-i+1?
        #       * For dest strip, is cell i identical to cell LEN-i+1?
        #       * For any 'Lookup' cells, does the strip have an odd number of cells and does the lookup index point to the center cell?

        return RewriteRule_Strip(
            Tuple(zip(lhs, rhs)),
            mask,
            weight,
            symmetries_explicit,
            symmetries_tail
        )
    finally
        pop!(inputs.op_stack_trace)
    end
end

"
Note that Source (lhs) patterns will return a Set as a List, to preserve input order for later processing.
The output element type is `Any` since it only loosely matches Source or Destination cell types.
"
function parse_markovjunior_rewrite_rule_md_side(inputs::MacroParserInputs,
                                                 loc, expr,
                                                 is_source::Bool
                                                )::Array{Any}
    with_parser_stacktrace(inputs,
                           is_source ?
                             "Left-hand side (multi-dimensional)" :
                             "Right-hand side (multi-dimensional)") do
        function lookup_char(c::Char, allow_wildcard::Bool)::Union{RewriteRuleCellSource, RewriteRuleCellDest{<:Tuple{Vararg{Int}}}}
            if haskey(CELL_CODE_BY_CHAR, c)
                CELL_CODE_BY_CHAR[c]
            elseif allow_wildcard && (c == '_')
                RewriteRuleCell_Wildcard()
            else
                raise_parse_error(
                    loc, inputs,
                    "Unsupported color value '", c, "'! Supported are ",
                    allow_wildcard ? "'_' and " : "",
                    "[ ", iter_join(keys(CELL_CODE_BY_CHAR), ", ")..., " ]"
                )
            end
        end
        function lookup_cell(cell_expr, pos)::Union{RewriteRuleCellSource, RewriteRuleCellDest{<:Tuple{Vararg{Int}}}}
            with_parser_stacktrace(inputs, "Cell $pos: `$cell_expr`") do
                return if @capture(cell_expr, {a_Symbol})
                    if is_source
                        raise_parse_error(loc, inputs,
                                          "Set syntax not allowed on the Source side -- use List syntax instead")
                    else
                        RewriteRuleCell_Set(string(a)...)
                    end
                elseif @capture(cell_expr, [a_Symbol])
                    RewriteRuleCell_List(Tuple(
                        lookup_char.(Iterators.map(identity, string(a)), Ref(false))
                    ))
                elseif cell_expr isa Symbol
                    cell_str = string(cell_expr)
                    if length(cell_str) != 1
                        raise_parse_error(loc, inputs,
                                          "Color char '", cell_str,
                                            "' should be one text character; got ", length(cell_str))
                    else
                        lookup_char(cell_str[1], true)
                    end
                elseif @capture(cell_expr, [a__])
                    if any(!isa(i, Int) for i in a)
                        raise_parse_error(loc, inputs,
                                          "Lookup index should be made of integers; got ", typeof.(a))
                    else
                        RewriteRuleCell_Lookup(Tuple(a))
                    end
                else
                    raise_parse_error(loc, inputs, "Unsupported syntax! `", cell_expr, "`")
                end
            end
        end

        # There are three kinds of multi-dimensional array literals:
        #  hcat (1D array) -- expr args are just the elements.
        #  vcat (2D array) -- expr args are each (head=:row, args=[ (values in this row....) ])
        #  ncat (3+D array) -- expr args are recursively (head=:ncat, args=[ dim, slice_exprs... ])
        #                        until getting to a 2D slice (head=:nrow, args=[ 1, row_exprs... ])
        #                        and each row is (head=:row, args=[ elements... ])
        as_array::Array{Any} = collect(Any, if Base.isexpr(expr, :hcat)
            # This snippet would normally generate a Vector,
            #    so make it a 1-row Matrix with 'permutedims()'.
            permutedims(lookup_cell.(expr.args, 1:length(expr.args)))
        elseif Base.isexpr(expr, :vcat)
            rows = (row_e -> row_e.args).(expr.args)
            raw_matrix = stack(rows, dims=1)
            lookup_cell.(raw_matrix, (c -> c.I).(CartesianIndices(raw_matrix)))
        elseif Base.isexpr(expr, :ncat)
            # Recursively stack rows, matrices, etc.
            function stack_arrays(slice_expr, pos::Tuple{Vararg{Int}})
                if Base.isexpr(slice_expr, :row)
                    # Note that these position tuples will treat X as the first axis,
                    #    which is fine as they're only used for logging errors to the user.
                    poses = (i -> (i, pos...)).(1:length(slice_expr.args))
                    return lookup_cell.(slice_expr.args, poses)
                elseif Base.isexpr(slice_expr, :nrow)
                    if isempty(slice_expr.args) || (slice_expr.args[1] != 1)
                        raise_parse_error(loc, inputs,
                            "Unexpected syntax format: :ncat/:nrow has ",
                            if isempty(slice_expr.args)
                                "no arguments"
                            else
                                "as its first argument `$(slice_expr.args[1])`"
                            end
                        )
                    end
                    elements = @view slice_expr.args[2:end]
                    poses = (i -> (i, pos...)).(1:length(elements))
                    return stack(stack_arrays.(elements, poses), dims=1)
                elseif Base.isexpr(slice_expr, :ncat)
                    if isempty(slice_expr.args) || !isa(slice_expr.args[1], Int)
                        raise_parse_error(loc, inputs,
                            "Unexpected syntax format: :ncat has ",
                            if isempty(slice_expr.args)
                                "no arguments"
                            else
                                "as its first argument `$(slice_expr.args[1])`"
                            end
                        )
                    end
                    elements = @view slice_expr.args[2:end]
                    poses = (i -> (i, pos...)).(1:length(elements))
                    return stack(stack_arrays.(elements, poses))
                else
                    raise_parse_error(loc, inputs,
                      "Unexpected syntax within :ncat, at slice ", pos, " -- ",
                      isa(slice_expr, Expr) ?
                        ":$(slice_expr.head)" :
                        "$(typeof(slice_expr))($slice_expr)"
                    )
                end
            end
            stack_arrays(expr, ())
        else
            error("Unexpected mode: ", expr)
        end)

        # Swap the first two axes, so that column (X) is the first instead of the second.
        # Note that 'list' rewrite cells do not need to be updated,
        #    as they are written assuming X as the first axis.
        # Also note that 1D arrays here will actually be a 1xN matrix, which still should do this
        as_array = permutedims(as_array, (2, 1, (3:ndims(as_array))...))

        return as_array
    end
end
function parse_markovjunior_rewrite_rule_md_symmetry(inputs::MacroParserInputs,
                                                     loc, exprs)::RewriteRule_MD_Symmetry_Definition
    # If no symmetry is given at all, then support all symmetries.
    if isnothing(exprs)
        return RewriteRule_MD_Symmetry_Definition()
    end

    with_parser_stacktrace(inputs, "symmetry statement") do
        s_permutations = Vector{Pair{Matrix{Int}, RewriteRule_TailSymmetry}}()
        s_chiralities = Vector{Set{Int}}()
        for expr in exprs
            # Look for a chirality statement.
            if @capture(expr, {chirality__})
                with_parser_stacktrace(inputs, "chirality statement $(length(s_chiralities) + 1)") do
                    new_axes = parse_markovjunior_axis_name.(Ref(loc), Ref(inputs), chirality)
                    unique_new_axes = Set(new_axes)
                    # I considered making it an error if there are redundant axes in this set,
                    #    but decided to prioritize usability over correctness.
                    if !isempty(unique_new_axes)
                        push!(s_chiralities, unique_new_axes)
                    end
                end
            # Look for a single-rule-axis statement.
            elseif @capture(expr, axis_Symbol[ values__ ])
                with_parser_stacktrace(inputs, "rule axis `$axis`") do
                    rule_axis = parse_markovjunior_axis_name(loc, inputs, axis)

                    # Parse out the normal permutations.
                    n_explicit_permutations = 0
                    grid_axes = map(values) do value_expr
                        if @capture(value_expr, iv_...)
                            Val(:tail) => parse_markovjunior_axis_signed(loc, inputs, iv)
                        else
                            (a, s) = parse_markovjunior_axis_signed(loc, inputs, value_expr)
                            n_explicit_permutations += 1
                            if s == 0
                                n_explicit_permutations += 1
                            end
                            (a, s)
                        end
                    end
                    permutation_matrix = fill(0,    1, n_explicit_permutations + 1)
                    permutation_matrix[1, 1] = rule_axis
                    n_written_permutations = 0
                    for element in grid_axes
                        if !isa(element, Pair{Val{:tail}})
                            (grid_axis, grid_sign) = element
                            if grid_sign in -1:0
                                n_written_permutations += 1
                                permutation_matrix[1, n_written_permutations + 1] = -grid_axis
                            end
                            if grid_sign in 0:1
                                n_written_permutations += 1
                                permutation_matrix[1, n_written_permutations + 1] = grid_axis
                            end
                        end
                    end

                    # Extract tail symmetries.
                    raw_tail_symmetry::NTuple{2, Optional{Int}} = (nothing, nothing)
                    for element in grid_axes
                        if element isa Pair{Val{:tail}}
                            (grid_axis, grid_sign) = element[2]
                            # Get the sign(s) covered by this tail symmetry,
                            #    and the index of each sign in the output.
                            signs::Tuple = if grid_sign == -1
                                tuple(1 => -1)
                            elseif grid_sign == 1
                                tuple(2 => 1)
                            elseif grid_sign == 0
                                tuple(1 => -1, 2 => 1)
                            else
                                error("Unhandled: ", element[2])
                            end

                            # Apply each sign.
                            for (s_i, grid_single_sign) in signs
                                if exists(raw_tail_symmetry[s_i])
                                    raise_parse_error(
                                        loc, inputs,
                                        "Conflicting tail symmetries given (second one is axis ",
                                          grid_axis, ")"
                                    )
                                else
                                    @set! raw_tail_symmetry[s_i] = grid_axis
                                end
                            end
                        end
                    end
                    # Simplify the data representation of tail symmetry.
                    tail_symmetry = if all(isnothing, raw_tail_symmetry)
                        nothing
                    elseif ==(raw_tail_symmetry...)
                        raw_tail_symmetry[1]
                    else
                        raw_tail_symmetry
                    end

                    push!(s_permutations, permutation_matrix => tail_symmetry)
                end
            # Look for a multi-rule-axis statement.
            elseif @capture(expr, rules_[ perms__ ]) && Base.isexpr(rules, :tuple)
                with_parser_stacktrace(inputs, "rule axes `$rules`") do
                    rule_axes::Vector{Int} = parse_markovjunior_axis_name.(Ref(loc), Ref(inputs), rules.args)
                    n_rule_axes = length(rule_axes)

                    # Each permutation given by the user can turn into one or more "real" permutations,
                    #    if they reference grid axes without giving a sign.
                    real_perms::Vector{Vector{Int}} = collect(Iterators.flatten(Iterators.map(perms) do perm_expr
                        with_parser_stacktrace(inputs, "Permutation $perm_expr") do
                            if Base.isexpr(perm_expr, :tuple)
                                perm_axes = parse_markovjunior_axis_signed.(Ref(loc), Ref(inputs), perm_expr.args)
                                if length(perm_axes) != n_rule_axes
                                    raise_parse_error(loc, inputs,
                                                      "Expected one element for each rule axis (",
                                                        n_rule_axes, "), got ",
                                                        length(perm_axes))
                                end

                                signed_output_perms = [ perm_axes ]
                                signed_perm_i = 1
                                while signed_perm_i <= length(signed_output_perms)
                                    signed_perm::Vector = signed_output_perms[signed_perm_i]
                                    for perm_axis_i in 1:length(signed_perm)
                                        (i_axis, i_sign) = signed_perm[perm_axis_i]
                                        if i_sign == 0
                                            signed_perm2 = copy(signed_perm)
                                            insert!(signed_output_perms, signed_perm_i + 1, signed_perm2)

                                            signed_perm[perm_axis_i] = (i_axis, -1)
                                            signed_perm2[perm_axis_i] = (i_axis, 1)

                                            signed_perm_i -= 1
                                            break
                                        end
                                    end

                                    signed_perm_i += 1
                                end

                                # Convert each permutation from tuple(axis, sign) to a signed integer.
                                return map(signed_output_perms) do perm
                                    return map(a -> a[1] * a[2], perm)
                                end
                            else
                                raise_parse_error(loc, inputs, "Expected a tuple of grid axes")
                            end
                        end
                    end))

                    # Build the properly-formatted permutation matrix.
                    n_permutations = length(real_perms)
                    permutation_matrix = fill(0, n_rule_axes, n_permutations + 1)
                    permutation_matrix[:, 1] = rule_axes
                    for perm_i in 1:n_permutations
                        permutation_matrix[:, perm_i + 1] = real_perms[perm_i]
                    end

                    # Note that multi-axis statements can't have a tail symmetry.
                    push!(s_permutations, permutation_matrix => nothing)
                end
            else
                raise_parse_error(
                    loc, inputs,
                    "Unexpected syntax: `", expr, "`! Expected one of: ",
                      "axis permutations (`x[ y, -z, +w... ]`), ",
                      "multi-axis permutations (`(x, z)[ (+x, +y) ]`), ",
                      "and chirality constraints (`{x, y, z}`)"
                )
            end
        end

        RewriteRule_MD_Symmetry_Definition(s_permutations, s_chiralities)
    end
end
function parse_markovjunior_rewrite_rule_md(inputs::MacroParserInputs, loc, expr)
    with_parser_stacktrace(inputs, "Block rewrite rule") do
        if !@capture(expr, lhsExpr_ => b_)
            raise_parse_error(loc, inputs, "Invalid format; expected `src => dest [modifiers]`")
        end

        # Strip out the modifiers.
        # NOTE: MacroTools has a bug where the | operator doesn't work.
        weightMulScalar = nothing
        weightDivScalar = nothing
        symmetryExprs = nothing
        maskExpr = nothing
        @capture(    b, c_ % maskExpr_ * weightMulScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ / weightDivScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ \ [ symmetryExprs__ ]) ||
            @capture(b, c_ % maskExpr_ * weightMulScalar_Real) ||
            @capture(b, c_ % maskExpr_ / weightDivScalar_Real) ||
            @capture(b, c_ * weightMulScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ / weightDivScalar_Real \[ symmetryExprs__ ]) ||
            @capture(b, c_ \[ symmetryExprs__ ]) ||
            @capture(b, c_ * weightMulScalar_Real) ||
            @capture(b, c_ / weightDivScalar_Real) ||
            @capture(b, c_ % maskExpr_) ||
            (c = b)
        rhsExpr = c

        # Parse the modifiers.
        weight = convert(Float32, if exists(weightMulScalar)
            @bp_check(isnothing(weightDivScalar))
            weightMulScalar
        elseif exists(weightDivScalar)
            1 / convert(Float64, weightDivScalar)
        else
            1
        end)
        mask = if @capture(maskExpr, maskAExpr_Real:maskBExpr_Real)
            convert.(Ref(Float32), (maskAExpr, maskBExpr))
        elseif maskExpr isa Real
            convert(Float32, maskExpr)
        elseif isnothing(maskExpr)
            nothing
        else
            raise_parse_error(loc, inputs, "Expected mask to be `%x` or `%(x:y)`; got `%$maskExpr`")
        end
        symmetry_def = parse_markovjunior_rewrite_rule_md_symmetry(inputs, loc, symmetryExprs)

        # Parse the rule.
        lhs = parse_markovjunior_rewrite_rule_md_side(inputs, loc, lhsExpr, true)
        rhs = parse_markovjunior_rewrite_rule_md_side(inputs, loc, rhsExpr, false)

        # Post-process and validate the rule.
        if size(lhs) != size(rhs)
            raise_parse_error(loc, inputs,
                              "Source block is ", iter_join(size(lhs), "x")...,
                                " while Dest block is ", iter_join(size(rhs), "x")...)
        end
        for ci::CartesianIndex in CartesianIndices(lhs)
            # If source cell was parsed as a list, remake it as a Set.
            # This reorders its elements, so if the destination is also a List then it needs to be reordered too.
            if lhs[ci] isa RewriteRuleCell_List
                if rhs[ci] isa RewriteRuleCell_List
                    if length(lhs[ci]) != length(rhs[ci])
                        raise_parse_error(loc, inputs,
                                          "Destination cell ", Tuple(ci), " is a list of ",
                                           length(rhs[ci]), " elements, but its Source cell has ",
                                           length(lhs[ci]), " elements instead")
                    else
                        lhs_to_rhs = Dict{UInt8, UInt8}(
                            L => R for (L, R) in zip(lhs[ci], rhs[ci])
                        )
                        lhs[ci] = RewriteRuleCell_Set(lhs[ci]...)
                        rhs[ci] = RewriteRuleCell_List(Tuple(lhs_to_rhs[L] for L in lhs[ci]))
                    end
                else
                    lhs[ci] = RewriteRuleCell_Set(lhs[ci]...)
                end
            end

            # If dest cell is a list, the source cell must be too (at this point actually a Set).
            if (rhs[ci] isa RewriteRuleCell_List) && !isa(lhs[ci], RewriteRuleCell_Set)
                raise_parse_error(loc, inputs,
                                  "Destination Cell ", Tuple(ci),
                                   " is a list, but its Source cell is not!")
            end

            # If dest cell references a source cell, that reference must be valid.
            if rhs[ci] isa RewriteRuleCell_Lookup{<:Tuple{Vararg{Int}}}
                source_idx = rhs[ci].source_idx
                if (length(source_idx) != ndims(lhs)) || any(!in(source_idx[n], 1:size(lhs, n))   for n in ndims(lhs))
                    raise_parse_error(loc, inputs,
                                      "Destination cell ", Tuple(ci),
                                        " references invalid source cell at ", source_idx)
                end
            end
        end

        # Make sure symmetry data is well-defined.
        #   * Rule axes aren't mentioned more than once.
        mentioned_rule_axes = fill(false, ndims(lhs))
        for (choice_set, tail_symmetry) in symmetry_def.grid_axis_choices
            for rule_axis in @view choice_set[:, 1]
                if mentioned_rule_axes[rule_axis]
                    raise_parse_error(loc, inputs,
                                      "More than one symmetry expression covers rule axis ",
                                      if rule_axis < 5
                                        ('x', 'y', 'z', 'w')[rule_axis]
                                      else
                                        rule_axis
                                      end)
                else
                    mentioned_rule_axes[rule_axis] = true
                end
            end
        end

        #TODO: Find symmetric axes (watch out for Lookup cells that reference one side!) and eliminate redundant explicit+tail-symmetrt options

        return RewriteRule_MD(
            collect(
                RewriteCell_MD{ndims(lhs)},
                Iterators.map(tuple, lhs, rhs)
            ),
            mask, weight,
            symmetry_def
        )
    end
end

function parse_markovjunior_rewrite_rules(inputs::MacroParserInputs, loc, expr
                                         )::Tuple{AbstractMarkovRewritePriority, Tuple{Vararg{Union{RewriteRule_Strip, RewriteRule_MD}}}}
    parse_rule_expr(inner_loc, line) = if @capture(line, a_ => b_) && any(Base.isexpr(a, n) for n in (:hcat, :vcat, :ncat)) 
        parse_markovjunior_rewrite_rule_md(inputs, inner_loc, line)
    else
        parse_markovjunior_rewrite_rule_strip(inputs, inner_loc, line)
    end
    if Base.isexpr(expr, :block)
        outputs = Vector{Union{RewriteRule_Strip, RewriteRule_MD}}()
        priority = Ref{Optional{AbstractMarkovRewritePriority}}(nothing)
        parse_markovjunior_block(expr.args) do inner_loc, line
            if @capture(line, PRIORITIZE(args__))
                with_parser_stacktrace(inputs, "PRIORITIZE(...)") do
                    if isempty(args)
                        raise_parse_error(loc, inputs, "Didn't provide a name!")
                    elseif exists(priority[])
                        raise_parse_error(loc, inputs, "Provided more than once!")
                    elseif !isa(args[1], Symbol)
                        raise_parse_error(loc, inputs, "The first argument isn't a name but a ", typeof(args[1]))
                    else
                        priority[] = parse_markovjunior_rewrite_priority(
                            Val(args[1]),
                            @view(args[2:end]),
                            inputs
                        )
                    end
                end
            else
                push!(outputs, parse_rule_expr(inner_loc, line))
            end
        end
        return (something(priority[], MarkovRewritePriority_Everything()),
                Tuple(outputs))
    else
        (MarkovRewritePriority_Everything(), tuple(parse_rule_expr(loc, expr)))
    end
end
function parse_markovjunior_op(::Val{Symbol("@rewrite")},
                               inputs::MacroParserInputs,
                               loc::Optional{LineNumberNode},
                               expr_args, original_line)
    args = Vector{Any}()
    parse_markovjunior_block(expr_args) do inner_loc, line
        push!(args, line)
    end

    if length(args) == 3
        with_parsed_markovjunior_bias_statement(inputs, loc, args[3]) do biases
            return MarkovOpRewrite(
                parse_markovjunior_rewrite_rules(inputs, loc, args[2])...,
                parse_markovjunior_threshold(inputs, loc, args[1]),
                Tuple(biases)
            )
        end
    elseif length(args) == 2
        if check_markovjunior_threshold_appearance(args[1])
            return MarkovOpRewrite(
                parse_markovjunior_rewrite_rules(inputs, loc, args[2])...,
                parse_markovjunior_threshold(inputs, loc, args[1]),
                ()
            )
        else
            with_parsed_markovjunior_bias_statement(inputs, loc, args[2]) do biases
                return MarkovOpRewrite(
                    parse_markovjunior_rewrite_rules(inputs, loc, args[1])...,
                    nothing,
                    Tuple(biases)
                )
            end
        end
    elseif length(args) == 1
        return MarkovOpRewrite(
            parse_markovjunior_rewrite_rules(inputs, loc, args[1])...,
            nothing,
            ()
        )
    elseif length(args) > 3
        raise_parse_error(nothing, inputs, "Should have 3 or fewer sections (threshold>rules>biases); got: ", args)
    elseif length(args) < 1
        raise_parse_error(nothing, inputs, "Statement is empty; needs to at least have one rewrite rule!")
    end
end


##############################
#  Rewrite Op: priority types

struct MarkovRewritePriority_Everything <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Everything,
                                          op::MarkovOpRewrite, state::MarkovOpRewrite_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    return weighted_random_array_element(
        (d.sum for d in state.weight_data_buffer_per_rule),
        state.weight_data_buffer.sum,
        rand(rng, Float32)
    )
end
function parse_markovjunior_rewrite_priority(::Val{:everything}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_parse_error(nothing, inputs,
                       "`everything` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Everything()
end
dsl_string(::MarkovRewritePriority_Everything) = "everything"

struct MarkovRewritePriority_Fair <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Fair,
                                          op::MarkovOpRewrite, state::MarkovOpRewrite_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    n_options = length(op.rules)
    chosen_i = rand(rng, 1:n_options)
    while isempty(rewrite_rule_option_indices(state, chosen_i)) # This rule has no matches?
        chosen_i = rand(rng, 1:n_options)
    end
    return chosen_i
end
function parse_markovjunior_rewrite_priority(::Val{:fair}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_parse_error(nothing, inputs,
                       "`fair` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Fair()
end
dsl_string(::MarkovRewritePriority_Fair) = "fair"

struct MarkovRewritePriority_Earliest <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Earliest,
                                          op::MarkovOpRewrite, state::MarkovOpRewrite_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    for i in 1:length(op.rules)
        if !isempty(rewrite_rule_option_indices(state, i))
            return i
        end
    end
    error("No rules have options???\n  Options = ", state.weighted_options_buffer,
             "\n  Idcs = ", state.weighted_options_buffer_first_indices)
end
function parse_markovjunior_rewrite_priority(::Val{:earliest}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_parse_error(nothing, inputs,
                       "`earliest` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Earliest()
end
dsl_string(::MarkovRewritePriority_Earliest) = "earliest"

struct MarkovRewritePriority_Latest <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Latest,
                                          op::MarkovOpRewrite, state::MarkovOpRewrite_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    for i in length(op.rules):-1:1
        if !isempty(rewrite_rule_option_indices(state, i))
            return i
        end
    end
    error("No rules have options???\n  Options = ", state.weighted_options_buffer,
             "\n  Idcs = ", state.weighted_options_buffer_first_indices)
end
function parse_markovjunior_rewrite_priority(::Val{:latest}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_parse_error(nothing, inputs,
                       "`latest` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Latest()
end
dsl_string(::MarkovRewritePriority_Latest) = "latest"

struct MarkovRewritePriority_Common <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Common,
                                          op::MarkovOpRewrite, state::MarkovOpRewrite_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    largest_i = 0
    largest_count::Float32 = -Inf32 # Float because counts are weighted
    for (i, rule) in zip(1:length(op.rules), op.rules)
        n_options = length(rewrite_rule_option_indices(state, i))
        next_count = n_options * rule.weight
        if next_count > largest_count
            largest_count = next_count
            largest_i = i
        end
    end
    return largest_i
end
function parse_markovjunior_rewrite_priority(::Val{:common}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_parse_error(nothing, inputs,
                       "`common` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Common()
end
dsl_string(::MarkovRewritePriority_Common) = "common"

struct MarkovRewritePriority_Rare <: AbstractMarkovRewritePriority end
function pick_rule_using_rewrite_priority(::MarkovRewritePriority_Rare,
                                          op::MarkovOpRewrite, state::MarkovOpRewrite_State,
                                          rng::PRNG, context::MarkovOpContext)::Int
    smallest_i = 0
    smallest_count::Float32 = Inf32 # Float because counts are weighted
    for (i, rule) in zip(1:length(op.rules), op.rules)
        n_options = length(rewrite_rule_option_indices(state, i))
        next_count = n_options * rule.weight
        if next_count < smallest_count
            smallest_count = next_count
            smallest_i = i
        end
    end
    return smallest_i
end
function parse_markovjunior_rewrite_priority(::Val{:rare}, expr_args, inputs::MacroParserInputs)
    if !isempty(expr_args)
        raise_parse_error(nothing, inputs,
                       "`rare` priority should have no arguments! Got ", length(expr_args))
    end
    return MarkovRewritePriority_Rare()
end
dsl_string(::MarkovRewritePriority_Rare) = "rare"