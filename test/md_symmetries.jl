# Wrap the test in a lambda to avoid name collisions.
(() -> begin

# Test find_all_md_symmetries()
function test_all_md_symmetries(n_rule_dims::Int, n_grid_dims::Int,
                                def::MJ.RewriteRule_MD_Symmetry_Definition,
                                expected::AbstractMatrix{MJ.GridDir}
                               )
    actual::AbstractMatrix{MJ.GridDir} = MJ.find_all_md_symmetries(def, n_rule_dims, n_grid_dims)
    n_entries = max(size(expected, 2), size(actual, 2))
    n_dims = max(size(expected, 1), size(actual, 1))

    function print_comparison()
        println("-------------------------------------------")
        println("Comparison |        Expected        Actual")
        println()
        for i_entry in 1:n_entries
            for i_dim in 1:n_dims
                val_expected = if (i_entry <= size(expected, 2)) && (i_dim <= size(expected, 1))
                    expected[i_dim, i_entry]
                else
                    nothing
                end
                val_actual = if (i_entry <= size(actual, 2)) && (i_dim <= size(actual, 1))
                    actual[i_dim, i_entry]
                else
                    nothing
                end

                print("   ")
                if i_dim == (n_dims÷2)
                    (i_entry < 10) && print(' ')
                    print(i_entry, ":")
                else
                    print("   ")
                end
                print("                  ")

                function print_val(v::Optional{MJ.GridDir})
                    if isnothing(v)
                        print(" ! ")
                    else
                        if v.axis < 10
                            print(' ')
                        end
                        print((v.sign > 0) ? '+' : '-')
                        if v.axis < 5
                            print(('x', 'y', 'z', 'w')[v.axis])
                        else
                            print(v.axis)
                        end
                    end
                end
                print_val(val_expected)
                print("          ")
                print_val(val_actual)

                println()
            end
            println("   ----    ---")
        end
        println("-------------------------------------------")
    end

    if size(actual) != size(expected)
        print_comparison()
        @bp_check(false,
            "\nERROR: Size mismatch! Expected ", size(expected), ", got ", size(actual),
            "\nSee printout above.\n"
        )
    end

    failed_msg = ""
    for i_entry in 1:n_entries
        for i_dim in 1:n_dims
            if actual[i_dim, i_entry] != expected[i_dim, i_entry]
                failed_msg *= string(
                    "Mismatch at entry ", i_entry, ", dimension ", i_dim, "!\n"
                )
            end
        end
    end
    if !isempty(failed_msg)
        print_comparison()
        @bp_check(false,
            "ERRORs: ", failed_msg, "\nSee printout above!"
        )
    end
end
test_all_md_symmetries(2, 2, MJ.RewriteRule_MD_Symmetry_Definition(
    [
     # R-Axis | Permutations...
        [ 1     2
          2     -1
        ] => nothing # <-- Tail symmetry
    ],
    # Chirality groups:
    Vector{Set{Int}}()
), reshape([  MJ.GridDir(2, 1), MJ.GridDir(1, -1) ], :, 1))
test_all_md_symmetries(2, 2, MJ.RewriteRule_MD_Symmetry_Definition(
    [
     # R-Axis | Permutations...
        [ 1     2    1
          2     -1   2
        ] => nothing # <-- Tail symmetry
    ],
    # Chirality groups:
    Vector{Set{Int}}()
    ), [ 
        MJ.GridDir(2, 1)   MJ.GridDir(1, 1)
        MJ.GridDir(1, -1)  MJ.GridDir(2, 1)
    ]
)
test_all_md_symmetries(2, 2, MJ.RewriteRule_MD_Symmetry_Definition(
    [
        # No explicit symmetries, but X has tail symmetry 2... and Y has +1...
        reshape([ 1 ], :, 1) => 2,
        reshape([ 2 ], :, 1) => (nothing, 1)
    ],
    Vector{Set{Int}}()
    ), [  MJ.GridDir(2, -1)   MJ.GridDir(2, 1)
        MJ.GridDir(1, 1)    MJ.GridDir(1, 1)
    ]
)
test_all_md_symmetries(1, 4, MJ.RewriteRule_MD_Symmetry_Definition(
    [
        reshape([ 1 ], :, 1) => (2, 4)
    ],
    Vector{Set{Int}}()
    ), reshape(
        [ MJ.GridDir(2, -1), MJ.GridDir(3, -1), MJ.GridDir(4, -1), MJ.GridDir(4, 1) ],
        1, :
    )
)
test_all_md_symmetries(4, 6, MJ.RewriteRule_MD_Symmetry_Definition(
    [
        [
            2      1    1    1
            4      -2   3    6
        ] => nothing,
        reshape([ 1 ], :, 1) => 5
    ],
    Vector{Set{Int}}()
    ), permutedims([ # Transposed so we can write each symmetry as a row
        # Fix 1 == {-5}
            # Fix 2,4 == {+1, -2}
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(3, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(3, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(6, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(6, +1)    MJ.GridDir(2, -1)
            # Fix 2,4 == {+1, +3}
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(2, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(2, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(6, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(6, +1)    MJ.GridDir(3, +1)
            # Fix 2,4 == {+1, +6}
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(2, -1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(2, +1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(3, -1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(3, +1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(6, +1)
        # Fix 1 == {+5}
            # Fix 2,4 == {+1, -2}
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(3, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(3, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(6, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(6, +1)    MJ.GridDir(2, -1)
            # Fix 2,4 == {+1, +3}
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(2, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(2, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(6, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(6, +1)    MJ.GridDir(3, +1)
            # Fix 2,4 == {+1, +6}
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(2, -1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(2, +1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(3, -1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(3, +1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(6, +1)
                MJ.GridDir(5, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(6, +1)
        # Fix 1 == {-6}
            # Fix 2,4 == {+1, -2}
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(3, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(3, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(5, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(5, +1)    MJ.GridDir(2, -1)
            # Fix 2,4 == {+1, +3}
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(2, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(2, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(5, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, -1)    MJ.GridDir(1, +1)    MJ.GridDir(5, +1)    MJ.GridDir(3, +1)
            # Fix 2,4 == {+1, +6}
                #    !  1 and 4 are both trying to be {6}  !
        # Fix 1 == {+6}
            # Fix 2,4 == {+1, -2}
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(3, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(3, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(5, -1)    MJ.GridDir(2, -1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(5, +1)    MJ.GridDir(2, -1)
            # Fix 2,4 == {+1, +3}
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(2, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(2, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(4, +1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(5, -1)    MJ.GridDir(3, +1)
                MJ.GridDir(6, +1)    MJ.GridDir(1, +1)    MJ.GridDir(5, +1)    MJ.GridDir(3, +1)
            # Fix 2,4 == {+1, +6}
                #    !  1 and 4 are both trying to be {6}  !
    ])
)
test_all_md_symmetries(2, 2, MJ.RewriteRule_MD_Symmetry_Definition(
    Pair{Matrix{Int}, MJ.RewriteRule_TailSymmetry}[ ],
    Set{Int}[ ]
), permutedims([
    MJ.GridDir(1, -1)  MJ.GridDir(2, -1)
    MJ.GridDir(1, -1)  MJ.GridDir(2, 1)
    MJ.GridDir(1, 1)   MJ.GridDir(2, -1)
    MJ.GridDir(1, 1)   MJ.GridDir(2, 1)
    MJ.GridDir(2, -1)  MJ.GridDir(1, -1)
    MJ.GridDir(2, -1)  MJ.GridDir(1, 1)
    MJ.GridDir(2, 1)  MJ.GridDir(1, -1)
    MJ.GridDir(2, 1)  MJ.GridDir(1, 1)
]))

# Test RewriteRule_MD_Orientations.
ori_ignored_fields(n_grid) = tuple(Vector{MJ.CellGrid{n_grid}}(), MJ.MarkovAllocatorHeap())
NULL_GRID_TO_RULE = MJ.GridDir(-1, -1, true)
# Note that your expected 'rule_permutation_cells_buffer' field can be left empty; we test its size internally.
# We also don't test anything about the allocator.
# orientations(::RewriteRule_MD_Orientations) is pretty bare-bones, so it's not tested either.
function test_orientations(src_array::Array{<:NTuple{2, Any}, NRule},
                           symmetry::MJ.RewriteRule_MD_Symmetry_Definition,
                           n_grid::Integer,
                           expected::MJ.RewriteRule_MD_Orientations) where {NRule}
    # If debug logging is enabled, print a clear separation from the previous test.
    if MJ.log_md_symmetry_logic()
        println(stderr, "\n\n\n===== test_orientations() =====")
    end

    sanitized_src_array = Array{MJ.RewriteCell_MD{NRule}, NRule}(undef, size(src_array))
    for i in eachindex(src_array)
        sanitize(cell) = if cell isa Integer
            convert(UInt8, cell)
        else
            cell
        end
        sanitized_src_array[i] = sanitize.(src_array[i]::Tuple)
    end

    actual = MJ.RewriteRule_MD_Orientations(sanitized_src_array, symmetry, Val(n_grid), MJ.MarkovAllocatorHeap())
    @bp_check(actual.rule_to_grid == expected.rule_to_grid,
        "Error getting MD-orientations! Rule-to-grid incorrect.\nRule: ",
          sprint(show, sanitized_src_array),
          "\nSymmetry: ", sprint(show, symmetry),
          "\nExpected: ", sprint(show, expected.rule_to_grid),
          "\nActual: ", sprint(show, actual.rule_to_grid)
    )
    @bp_check(actual.grid_to_rule == expected.grid_to_rule,
        "Error getting MD-orientations! Grid-to-rule incorrect.\nRule: ",
          sprint(show, sanitized_src_array),
          "\nSymmetry: ", sprint(show, symmetry),
          "\nExpected: ", sprint(show, expected.grid_to_rule),
          "\nActual: ", sprint(show, actual.grid_to_rule)
    )
    @bp_check(length(actual.rule_permutations) == length(expected.rule_permutations),
        "Error getting MD-orientations! `rule_permutations` has a mismatch: ",
          "expected ", length(expected.rule_permutations), " but got ", length(actual.rule_permutations)
    )
    @bp_check(actual.rule_permutations == expected.rule_permutations,
        "Error getting MD-orientations! `rule_permutations` has a mismatch.\nRule: ",
          sprint(show, sanitized_src_array),
          "\nSymmetry: ", sprint(show, symmetry),
          "\nExpected:\n\t", join(sprint.(Ref(show), expected.rule_permutations), "\n\t"),
          "\nActual:\n\t", join(sprint.(Ref(show), actual.rule_permutations), "\n\t")
    )
    @bp_check(size.(actual.rule_permutation_cells_buffer) == size.(expected.rule_permutations),
        "Error getting MD-orientations! `rule_permutation_cells_buffer` elements have the wrong sizes.",
          "\nRule: ", sprint(show, sanitized_src_array),
          "\nSymmetry: ", sprint(show, symmetry),
          "\nExpected:\n\t", size.(expected.rule_permutation),
          "\nActual:\n\t", size.(actual.rule_permutation_cells_buffer)
    )
end
test_orientations(
    [  # (remember, first axis is Row)
        (2,4)  (4,8)
        (3, MJ.RewriteRuleCell_Lookup((1, 2)))  (5,10)
    ],
    MJ.RewriteRule_MD_Symmetry_Definition(
        [
           # Rule-Axis | Permutations...
            [ 1          2
              2         -1
            ] => nothing # <-- Tail symmetry
        ],
        # Chirality groups:
        Vector{Set{Int}}()
    ),
    2,
    # Expected:
    MJ.RewriteRule_MD_Orientations(
        # Expect only one permutation: [ 4 2 ;; 5 3 ]
        1,
        permutedims([ MJ.GridDir(2, 1)    MJ.GridDir(1, -1) ]),
        permutedims([ MJ.GridDir(2, -1)   MJ.GridDir(1, 1) ]),

        Array{MJ.RewriteCell_MD{2}, 2}[
            collect(MJ.RewriteCell_MD{2}, permutedims([
                (0x04, 0x08)  (0x02, 0x04)
                (0x05, 0x0A)   (0x03, MJ.RewriteRuleCell_Lookup((1, 1)))
            ]))
        ],
        ori_ignored_fields(2)...
    )
)
test_orientations(
    [ (1,2), (3, MJ.RewriteRuleCell_Lookup((3, ))), (5,6) ],
    MJ.RewriteRule_MD_Symmetry_Definition(
        # The rule can only be applied in reverse on the X axis, but any way on axes Z and up.
        [
            [ 1   -1 ] => 3
        ],
        Vector{Set{Int}}()
    ),
    4,
    MJ.RewriteRule_MD_Orientations(
        # -X, -Z, +Z, -W, +W
        5,
        [ MJ.GridDir(1, -1)   MJ.GridDir(3, -1)   MJ.GridDir(3, 1)   MJ.GridDir(4, -1)   MJ.GridDir(4, 1) ],
        permutedims([ # Each row here becomes a column (orientation enry) in the test data.
            MJ.GridDir(1, -1)  NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  NULL_GRID_TO_RULE
            NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  MJ.GridDir(1, -1)  NULL_GRID_TO_RULE
            NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  MJ.GridDir(1, 1)  NULL_GRID_TO_RULE
            NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  MJ.GridDir(1, -1)
            NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  NULL_GRID_TO_RULE  MJ.GridDir(1, 1)
        ]),

        Array{MJ.RewriteCell_MD{4}, 4}[
            collect.(Ref(MJ.RewriteCell_MD{4}), [
                # -X
                [ UInt8.((5,6)) ; (UInt8(3), MJ.RewriteRuleCell_Lookup((1, 1, 1, 1))) ; UInt8.((1,2)) ;;;; ],
                # -Z, +Z
                [ UInt8.((5,6)) ;;; (UInt8(3), MJ.RewriteRuleCell_Lookup((1, 1, 1, 1))) ;;; UInt8.((1,2)) ;;;; ],
                [ UInt8.((1,2)) ;;; (UInt8(3), MJ.RewriteRuleCell_Lookup((1, 1, 3, 1))) ;;; UInt8.((5,6)) ;;;; ],
                # -W, +W
                [ UInt8.((5,6)) ;;;; (UInt8(3), MJ.RewriteRuleCell_Lookup((1, 1, 1, 1))) ;;;; UInt8.((1,2)) ],
                [ UInt8.((1,2)) ;;;; (UInt8(3), MJ.RewriteRuleCell_Lookup((1, 1, 1, 3))) ;;;; UInt8.((5,6)) ]
            ])...
        ],
        ori_ignored_fields(4)...
    )
)
test_orientations(
    [
        (1,2)    (MJ.RewriteRuleCell_Wildcard(), 4)    (5, MJ.RewriteRuleCell_Wildcard())
        (MJ.RewriteRuleCell_Set(0x9, 0xa), MJ.RewriteRuleCell_List((0x7, 0x8)))  (11,12)  (13,14)
    ],
    MJ.RewriteRule_MD_Symmetry_Definition(
        [
            [ 1     1  -2  3
              2     2   1  -1
            ] => nothing
        ],
        Vector{Set{Int}}()
    ),

    3,
    MJ.RewriteRule_MD_Orientations(
        3,
        [ MJ.GridDir(1, 1)    MJ.GridDir(2, -1)   MJ.GridDir(3, 1)
          MJ.GridDir(2, 1)    MJ.GridDir(1, 1)    MJ.GridDir(1, -1)
        ],
        [
            MJ.GridDir(1, 1)   MJ.GridDir(2, 1)    MJ.GridDir(2, -1)
            MJ.GridDir(2, 1)   MJ.GridDir(1, -1)   NULL_GRID_TO_RULE
            NULL_GRID_TO_RULE  NULL_GRID_TO_RULE   MJ.GridDir(1, 1)
        ],

        Array{MJ.RewriteCell_MD{3}, 3}[
            collect.(Ref(MJ.RewriteCell_MD{3}), [
                [
                    (0x1, 0x2)    (MJ.RewriteRuleCell_Wildcard(), 0x4)    (0x5, MJ.RewriteRuleCell_Wildcard())
                    (MJ.RewriteRuleCell_Set(0x9, 0xa), MJ.RewriteRuleCell_List((0x7, 0x8)))  (0xb, 0xc)  (0xd, 0xe)
                    ;;;
                ],
                [
                    (MJ.RewriteRuleCell_Set(0x9, 0xa), MJ.RewriteRuleCell_List((0x7, 0x8)))   (0x1, 0x2)
                    (0xb, 0xc)   (MJ.RewriteRuleCell_Wildcard(), 0x4)
                    (0xd, 0xe)   (0x5, MJ.RewriteRuleCell_Wildcard())
                    ;;;
                ],
                [
                    (0x5, MJ.RewriteRuleCell_Wildcard())
                    (MJ.RewriteRuleCell_Wildcard(), 0x4)
                    (0x1, 0x2)

                    ;;;

                    (0xd, 0xe)
                    (0xb, 0xc)
                    (MJ.RewriteRuleCell_Set(0x9, 0xa), MJ.RewriteRuleCell_List((0x7, 0x8)))
                ]
            ])...
        ],
        ori_ignored_fields(3)...
    )
)

end)()