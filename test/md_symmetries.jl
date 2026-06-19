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

end)()