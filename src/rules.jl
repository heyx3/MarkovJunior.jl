"
An input (matching some cell line in a grid), and an output (replacing that line).
Rules are one-dimensional, and can be applied in either direction along any axis.

Null entries in the input/output are wildcards --
    wild inputs don't care about that cell, and wild outputs leave that cell unchanged.

On construction you can specify input and output as strings, or arrays of char/int/Nothing.
"
struct CellRule
    length::Int32
    input::Vector{Optional{UInt8}}
    output::Vector{Optional{UInt8}}

    function CellRule(input::Union{AbstractVector, AbstractString, Symbol},
                      output::Union{AbstractVector, AbstractString, Symbol})
        if input isa Symbol || output isa Symbol
            return CellRule(input isa Symbol ? string(input) : input,
                            output isa Symbol ? string(output) : output)
        end

        if length(input) != length(output)
            error("Input and output rules must be the same size!")
        end
        converter = list -> collect(
            Optional{UInt8},
            Iterators.map(list) do i
                if i isa AbstractChar
                    c = convert(Char, i)
                    if c in keys(CELL_CODE_BY_CHAR)
                        CELL_CODE_BY_CHAR[convert(Char, i)]
                    elseif c == CELL_CHAR_INVALID
                        nothing
                    else
                        error("Unsupported rule char '$c', expected one of ",
                              collect(keys(CELL_CODE_BY_CHAR)))
                    end
                elseif i isa Integer
                    convert(UInt8, i)
                elseif isnothing(i)
                    nothing
                else
                    error("Unhandled type of rule cell: ", typeof(i))
                end
            end
        )
        return new(convert(Int32, length(input)),
                   converter(input), converter(output))
    end
end

dsl_string(rule_string::Vector{Optional{UInt8}}) = String(map(r -> isnothing(r) ? CELL_CHAR_INVALID : dsl_string(r), rule_string))
dsl_string(r::CellRule) = "@rule $(dsl_string(r.input)) => $(dsl_string(r.output))"

function rule_applies(grid::CellGrid{N},
                      rule::CellRule,
                      at::CellLine{N}
                     )::Bool where {N}
    @markovjunior_assert at.length == rule.length
    success = true
    for_each_cell_in_line(at) do i::Int32, cell::CellIdx{N}
        if exists(rule.input[i])
            success &= rule.input[i] == grid[cell]
        end
    end
    return success
end
function rule_execute!(grid::CellGrid{N},
                       rule::CellRule,
                       at::CellLine{N}
                      )::Nothing where {N}
    @markovjunior_assert at.length == rule.length
    for_each_cell_in_line(at) do i::Int32, cell::CellIdx{N}
        if exists(rule.output[i])
            grid[cell] = rule.output[i]
        end
    end
    return nothing
end

"
Finds all matches of the given rule in the given grid, sending them into your lambda.

This is mainly for initialization and debugging;
    you should use a `RuleCache` for better performance across multiple rule applications.
"
function find_rule_matches(process_rule, # CellLine{N} -> Nothing
                           grid::CellGrid{N},
                           rule::CellRule
                          )::Nothing where {N}
    for dir_idx in 1:grid_dir_count(N)
        dir = grid_dir_index(dir_idx)

        # If a rule is length-1, then different grid directions are redundant.
        (rule.length <= 1 && dir_idx > 1) && continue

        first_grid_idx::CellIdx{N} = one(CellIdx{N})
        last_grid_idx::CellIdx{N} = vsize(grid)
        if dir.dir < 0
            @set! first_grid_idx[dir.axis] += rule.length - 1
        else
            @set! last_grid_idx[dir.axis] -= rule.length - 1
        end

        for grid_idx::CellIdx{N} in first_grid_idx:last_grid_idx
            at = CellLine(grid_idx, dir, rule.length)
            if rule_applies(grid, rule, at)
                process_rule(at)
            end
        end
    end
    return nothing
end
"
Find matches of the given rules in the given grid, and send them into your lambda.
Rules are ordered by your own rules iterator,
    e.g. the first batch of outputs are all for your first rule.

This is mainly for initialization and debugging;
    you should use a `RuleCache` for better performance across multiple rule applications.
"
function find_rule_matches(process_rule, # (rule_idx::Int32, CellLine{N}) -> Nothing
                           grid::CellGrid{N},
                           rules::AbstractVector{CellRule}
                          )::Nothing where {N}
    for (i, rule) in enumerate(rules)
        find_rule_matches(at -> process_rule(convert(Int32, i), at), grid, rule)
    end
    return nothing
end



"
Stored lookup tables indicating which rules are currently applicable in which cells of some grid.
You must always call `update_cache!()` after updating the grid.
"
struct RuleCache{N}
    # For each rule, stores all the possible legal applications of it
    #    as a tuple of (first cell idx, direction).
    #
    # The set of legal applications is ordered so that we can deterministically iterate over it.
    legal_applications::Vector{OrderedSet{Tuple{CellIdx{N}, GridDir}}}
    # Counts how many cached legal rule applications affect each cell.
    count_per_cell::Array{Int32, N}

    rules::Vector{CellRule}


    function RuleCache(grid::CellGrid{N}, rules::AbstractVector{CellRule}) where {N}
        legal_applications = map(r -> OrderedSet{Tuple{CellIdx{N}, GridDir}}(), rules)
        count_per_cell = fill(zero(Int32), size(grid))
        find_rule_matches(grid, rules) do rule_idx::Int32, at::CellLine{N}
            push!(legal_applications[rule_idx], (at.start_cell, at.movement))
            for_each_cell_in_line(at) do i::Int32, cell::CellIdx{N}
                count_per_cell[cell] += 1
            end
        end

        return new{N}(legal_applications, count_per_cell, collect(rules))
    end
end


"
Updates the rule cache after the given rule application.
This must always be called immediately after any changes to the grid.
"
function update_cache!(rc::RuleCache{N}, grid::CellGrid{N},
                       executed_rule_idx::Int32,
                       executed_rule_at::CellLine{N}
                      )::Nothing where {N}
    @markovjunior_assert(
        contains(box_indices(grid), cell_line_aabb(executed_rule_at)),
        "Rule application passed to update_cache!() doesn't actually fit in the grid"
        #NOTE: some logic below does not clamp cell indices because of this assumption
    )
    executed_rule = rc.rules[executed_rule_idx]
    executed_axis = executed_rule_at.movement.axis

    executed_min_parallel = executed_rule_at.start_cell[executed_rule_at.movement.axis]
    executed_max_parallel = executed_min_parallel + (executed_rule_at.movement.dir * (executed_rule.length - 1))
    (executed_min_parallel, executed_max_parallel) = minmax(executed_min_parallel, executed_max_parallel)

    # Go through every possible rule application that involves one of the cells in this line,
    #    and re-evaluate whether it fits.

    function evaluate_rule_application(rule_idx::Int32, at_start::CellIdx{N}, at_dir::GridDir)
        rule = rc.rules[rule_idx]
        at = CellLine{N}(at_start, at_dir, rc.rules[rule_idx].length)
        application_key = (at_start, at_dir)

        used_to_be_legal::Bool = application_key in rc.legal_applications[rule_idx]
        is_legal::Bool = rule_applies(grid, rc.rules[rule_idx], at)

        # Became legal?
        if !used_to_be_legal && is_legal
            push!(rc.legal_applications[rule_idx], application_key)
            for_each_cell_in_line(at) do idx, cell::CellIdx{N}
                rc.count_per_cell[cell] += 1
            end
        # Became illegal?
        elseif used_to_be_legal && !is_legal
            delete!(rc.legal_applications[rule_idx], application_key)
            for_each_cell_in_line(at) do idx, cell::CellIdx{N}
                rc.count_per_cell[cell] -= 1
            end
        end
    end

    for (_rule_idx, rule) in enumerate(rc.rules)
        rule_idx = convert(Int32, _rule_idx)

        # First evaluate rules ahead and behind this cell line, intersecting it through either end.
        for parallel_sign in Int32.((-1, +1))
            # Only one 'sign' is needed in the case of a length-1 rule.
            (rule.length <= 1 && parallel_sign == 1) && break

            parallel_grid_size = size(grid, executed_axis)
            (parallel_min, parallel_max) = if parallel_sign < 0
                (max(executed_min_parallel,
                     rule.length),
                 min(executed_max_parallel + rule.length - one(Int32),
                     parallel_grid_size))
            else
                (max(executed_min_parallel - rule.length + one(Int32),
                     one(Int32)),
                 min(executed_max_parallel,
                     parallel_grid_size - rule.length + 1))
            end

            for parallel_pos in parallel_min:parallel_max
                rule_start = executed_rule_at.start_cell
                @set! rule_start[executed_axis] = parallel_pos

                # If the rule is length-1, then all grid directions are redundant
                #    and we standardize on the first axis.
                reported_axis = if rule.length <= 1
                    one(Int32)
                else
                    executed_axis
                end
                evaluate_rule_application(rule_idx, rule_start,
                                          GridDir(reported_axis, parallel_sign))
            end
        end

        # Next evaluate rules perpendicular to this cell line, intersecting it at exactly one cell.
        # This work is redundant in the case of a length-1 rule.
        (rule.length > 1) && for perp_axis_idx in one(Int32):convert(Int32, N)
            if perp_axis_idx == executed_axis
                continue
            end
            perp_grid_size = size(grid)[perp_axis_idx]
            perp_focus = executed_rule_at.start_cell[perp_axis_idx]

            # NOTE: no need to clamp parallel position, as we asserted this whole line is in the grid.
            for parallel_offset in zero(Int32):convert(Int32, executed_rule.length - 1)
                for perp_dir in Int32.((-1, +1))
                    (perp_start_min, perp_start_max) = if perp_dir > 0
                        (max(one(Int32),
                             perp_focus - rule.length + one(Int32)),
                         min(perp_grid_size - rule.length + one(Int32),
                             perp_focus))
                    else
                        (max(one(Int32) + rule.length - one(Int32),
                             perp_focus),
                         min(perp_grid_size,
                             perp_focus + rule.length - one(Int32)))
                    end

                    for perp_pos in perp_start_min:perp_start_max
                        rule_start = executed_rule_at.start_cell
                        @set! rule_start[executed_axis] += (parallel_offset * executed_rule_at.movement.dir)
                        @set! rule_start[perp_axis_idx] = perp_pos

                        rule_dir = GridDir(perp_axis_idx, perp_dir)
                        evaluate_rule_application(rule_idx, rule_start, rule_dir)
                    end
                end
            end
        end
    end
end


"Calculates the total number of possible rule applications"
n_cached_rule_applications(c::RuleCache) = Int32(sum(Iterators.map(length, c.legal_applications), init=0))


"An entry in the rule cache, representing a particular rule applied to particular cells"
struct RuleApplication{N}
    rule_idx::Int32
    at::CellLine{N}
end

"
Runs your lambda on all possible rule applications.
If you want to grab a specific one, `get_cached_rule_application()` is more efficient.
"
function for_each_cached_rule_application(lambda, c::RuleCache)
    for (rule_idx, cell_lines) in enumerate(c.legal_applications)
        for (start_pos, dir) in cell_lines
            lambda(RuleApplication(
                convert(Int32, rule_idx),
                CellLine(start_pos, dir, c.rules[rule_idx].length)
            ))
        end
    end
end

"
Gets a specific legal rule application from the cache, by its index.
If you want to iterate over all of them, `for_each_cached_rule_application()` is more efficient.
"
function get_cached_rule_application(c::RuleCache{N}, i::Int)::RuleApplication{N} where {N}
    @bp_check(i > 0, "Index must be at least 1; got ", i)

    relative_i = i
    rule_idx::Int = 1
    while (rule_idx <= length(c.rules)) && (relative_i > length(c.legal_applications[rule_idx]))
        relative_i -= length(c.legal_applications[rule_idx])
        rule_idx += 1
    end
    @markovjunior_assert(relative_i > 0)
    @bp_check(rule_idx <= length(c.rules),
              "Index ", i, " out of range 1:", n_cached_rule_applications(c))

    (cell_start, movement) = ordered_collection_get_idx(c.legal_applications[rule_idx], relative_i)
    return RuleApplication{N}(rule_idx, CellLine{N}(cell_start, movement, c.rules[rule_idx].length))
end

"
Gets the set of all possible applications, for the first rule which has at least one.
Also returns that rule, by index.

**The returned set is read-only!**

If no options are available, returns nothing.
"
function get_first_nonempty_cached_rule_application_set(
             c::RuleCache{N}
         )::Optional{Tuple{Int32, OrderedSet{Tuple{CellIdx{N}, GridDir}}}} where {N}
    for rule_idx in 1:length(c.rules)
        if !isempty(c.legal_applications[rule_idx])
            return (convert(Int32, rule_idx), c.legal_applications[rule_idx])
        end
    end
    return nothing
end