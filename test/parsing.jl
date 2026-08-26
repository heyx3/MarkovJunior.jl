# Wrap the test in a lambda to avoid name collisions.
(() -> begin

# Define helpers to compare two parsed algorithms field-by-field.
function test_compare(a::MJ.MarkovAlgorithm, b::MJ.MarkovAlgorithm, tab::String = "")
    println(tab, "Comparing two algorithms...")

    if length(a.sequence) != length(b.sequence)
        println(tab, "\tA has ", length(a.sequence), " elements while B has ", length(b.sequence), "!")
    else
        for (oa, ob, i) in zip(a.sequence, b.sequence, 1:length(a.sequence))
            println(tab, "\tOp ", i, ", ", typeof(oa), "(", typeof(ob), ")...")
            test_compare(oa, ob, "$tab\t\t")
        end
    end

    if (length(a.pragmas_chronological) != length(b.pragmas_chronological)) ||
       (length(a.pragmas_map) != length(b.pragmas_map))
    #begin
       println(tab, "\tPragmas don't line up! A has ",
                    length(a.pragmas_chronological), "/", length(a.pragmas_map),
                    " and B has ", length(b.pragmas_chronological), "/", length(b.pragmas_map))
    else
        for ((pan, pai), (pbn, pbi), i) in zip(a.pragmas_chronological, b.pragmas_chronological, 1:length(a.pragmas_chronological))
            if (pan, pai) != (pbn, pbi)
                println(tab, "\tPragma ", i, " is different! A is ",
                             (pan, pai), " vs B is ", (pbn, pbi))
            end
        end
        for (a_name, a_elements) in a.pragmas_map
            !haskey(b.pragmas_map, a_name) && continue # Would already have been logged above
            b_elements = b.pragmas_map[a_name]
            if length(a_elements) != length(b_elements)
                println(tab, "\tPragmas under ", a_name, " don't line up!",
                             "A has ", length(a_elements), "entries while B has ", length(b_elements),
                             ":\n", tab, "\t\tA=", a_elements,
                             "\n", tab, "\t\tB=", b_elements)
            else
                for (a_entry, b_entry, entry_i) in zip(a_elements, b_elements, 1:length(a_elements))
                    if a_entry != b_entry
                        println(tab, "\tPragma[", a_name, "][", entry_i, "] doesn't line up! ",
                                     "A has ", a_entry, " vs B has ", b_entry)
                    end
                end
            end
        end
    end

    # Note: 'add_ons' is a runtime parameter and not worth printing here.

    return nothing
end

test_compare(a::MJ.AbstractMarkovOp, b::MJ.AbstractMarkovOp, tab::String) = println(tab, "Mismatched/Unsupported types!")
function test_compare(a::MJ.MarkovOpSequence, b::MJ.MarkovOpSequence, tab::String)
    if a.threshold != b.threshold
        println(tab, "Mismatch of thresholds! A has ", a.threshold, " while B has ", b.thresold)
    end
    if length(a.ops) != length(b.ops)
        println(tab, "A has ", length(a.ops), " inner ops, but B has ", length(b.ops), "!")
    else
        for (oa, ob, i) in zip(a.ops, b.ops, 1:length(a.ops))
            println(tab, "Op ", i, ", ", MJ.dsl_string(oa), "   (", MJ.dsl_string(ob), ")")
            test_compare(oa, ob, "$tab\t")
        end
        if a.ops != b.ops
            println(tab, "Op list Mismatch!")
        end
    end
    if length(a.biases) != length(b.biases)
        println(tab, "A has ", length(a.biases), " biases, but B has ", length(b.biases), "!")
    else
        for (ba, bb, i) in zip(a.biases, b.biases, 1:length(a.biases))
            println(tab, "Bias ", i)
            test_compare(ba, bb, "$tab\t")
        end
    end
    if a != b
        println(tab, "Fails to match via == operator!")
    end
end
function test_compare(a::MJ.MarkovOpDrawBox, b::MJ.MarkovOpDrawBox, tab::String)
    for f in fieldnames(typeof(a))
        af = getfield(a, f)
        bf = getfield(b, f)
        if af != bf
            println(tab, "Mismatch of ", f, "! A has ", af, " while B has ", bf)
        end
    end
    if a != b
        println(tab, "Fails to match via == operator!")
    end
end
function test_compare(a::MJ.MarkovOpRewrite, b::MJ.MarkovOpRewrite, tab::String)
    if length(a.rules) != length(b.rules)
        println(tab, "A has ", length(a.rules), " rules while B has ", length(b.rules), "!")
    else
        for (ra, rb, i) in zip(a.rules, b.rules, 1:length(a.rules))
            println(tab, "Rule ", i, ", ", MJ.dsl_string(ra), "    (", MJ.dsl_string(rb), ")")
            if typeof(ra).name != typeof(rb).name
                println(tab, "\tA is ", typeof(ra).name, " while B is ", typeof(rb).name, "!")
                return nothing
            end

            cell_count_fn = (ra isa MJ.RewriteRule_Strip) ? length : size
            cell_vsize_fn = (ra isa MJ.RewriteRule_Strip) ? (t -> Vec(length(t))) : vsize
            if cell_count_fn(ra.cells) != cell_count_fn(rb.cells)
                println(tab, "\tA has ", iter_join(cell_count_fn(ra.cells), "x")...,
                             " cells, while B has ", iter_join(cell_count_fn(rb.cells), "x")...,
                             "!")
            else
                v_n_cells = cell_vsize_fn(ra.cells)
                for (ca, cb, _j) in zip(ra.cells, rb.cells, one(typeof(v_n_cells)):v_n_cells)
                    j = if ndims(_j) == 1
                        _j.x
                    else
                        _j
                    end
                    println(tab, "\tCell ", j, ", ", ca, "   (", cb, ")")
                    if ca != cb
                        println(tab, "\t\tMismatch at cell ", j, "!\n",
                                tab, "\t\t\tA=", ca, "\n", tab, "\t\t\tB=", cb)
                    end
                end
            end
            if ra.weight != rb.weight
                println(tab, "\tWeight mismatch! ", ra.weight, " vs ", rb.weight)
            end
            if ra.mask != rb.mask
                println(tab, "\tMask mismatch! ",
                        typeof(ra.mask), "(", ra.mask, ") vs ",
                        typeof(rb.mask), "(", rb.mask, ")")
            end
            if (ra isa MJ.RewriteRule_Strip)
                if length(ra.explicit_symmetries) != length(rb.explicit_symmetries)
                    println(tab, "\tA has ", length(ra.explicit_symmetries), " explicit symmetries ",
                            "while B has ", length(rb.explicit_symmetries))
                else
                    for (sa, sb, j) in zip(ra.explicit_symmetries, rb.explicit_symmetries, 1:length(ra.explicit_symmetries))
                        println(tab, "\tExplicit symmetry ", sa, "  (", sb, ")...")
                        if sa != sb
                            println(tab, "\t\tMismatch!")
                        end
                    end
                end
                if ra.tail_symmetry != rb.tail_symmetry
                    println(tab, "\tTail-symmetry mismatch! ",
                                ra.tail_symmetry,
                                " vs ", rb.tail_symmetry)
                end
            elseif (ra isa MJ.RewriteRule_MD)
                (sa, sb) = (ra.symmetry, rb.symmetry)
                if length(sa.grid_axis_choices) != length(sb.grid_axis_choices)
                    println(tab, "\tSymmetry mismatch! ",
                                 "A has ", length(sa.grid_axis_choices), " axis choice sets ",
                                 "while B has ", length(sb.grid_axis_choices))
                else for ((ma, ta), (mb, tb), mi) in zip(sa.grid_axis_choices, sb.grid_axis_choices, 1:length(sa.grid_axis_choices))
                    if @view(ma[:, 1]) != @view(mb[:, 1])
                        println(tab, "\tMismatch of grid_axis_choices[", mi, "]'s chosen rule axes!\n",
                                tab, "\t\tA=", ma[:, 1], "\n",
                                tab, "\t\tB=", mb[:, 1])
                    end
                    if @view(ma[:, 2:end]) != @view(mb[:, 2:end])
                        println(tab, "\tMismatch of grid_axis_choices[", mi, "]'s permutations!\n",
                                tab, "\t\tA: ", ma[:, 2:end], "\n",
                                tab, "\t\tB: ", mb[:, 2:end])
                    end
                    if ta != tb
                        println(tab, "\tMismatch of grid_axis_choices[", mi, "]'s tail symmetry!\n",
                                tab, "\t\tA=", ta, "   B=", tb)
                    end
                end end
                if length(sa.chiral_groups) != length(sb.chiral_groups)
                    println(tab, "\tSymmetry mismatch! ",
                                 "A has ", length(sa.chiral_groups), " chiral groups ",
                                 "while B  has ", length(sb.chiral_groups))
                else for (ca, cb, ci) in zip(sa.chiral_groups, sb.chiral_groups, 1:length(sa.chiral_groups))
                    println(tab, "\tSymmetry chiral group ", ci, ": A={", iter_join(ca, ", ")...,
                                 "}    B={", iter_join(cb, ", ")..., "}")
                    if ca != cb
                        println(tab, "\t\tMismatch! ")
                    end
                end end
            else
                error("Unhandled: ", typeof(ra))
            end

            if ra != rb
                println(tab, "\tFails to match via == operator!")
            end
        end
    end
    if a.threshold != b.threshold
        println(tab, "\tThreshold mismatch! ",
                typeof(a.threshold), "(", a.threshold, ") vs ",
                typeof(b.threshold), "(", b.threshold, ")")
    end
    if length(a.biases) != length(b.biases)
        println(tab, "A has ", length(a.biases), " biases, but B has ", length(b.biases), "!")
    else
        for (ba, bb, i) in zip(a.biases, b.biases, 1:length(a.biases))
            println(tab, "Bias ", i)
            test_compare(ba, bb, "$tab\t")
        end
        if a.biases != b.biases
            println(tab, "Bias tuple Mismatch!")
        end
    end
    if a != b
        println(tab, "Fails to match via == operator!")
    end
    return nothing
end

function test_compare(a::MJ.AbstractMarkovBias, b::MJ.AbstractMarkovBias, tab::String)
    if typeof(a) != typeof(b)
        println(tab, "Mismatched/Unsupported types! A is ", typeof(a), " while B is ", typeof(b))
    else
        println(tab, "A \"", MJ.dsl_string(a), "\"   vs B \"", MJ.dsl_string(b), "\"")
        for p_name in propertynames(a)
            p_a = getproperty(a, p_name)
            p_b = getproperty(b, p_name)
            if p_a != p_b
                println(tab, "Mismatched ", p_name, "! A has ", p_a, " while B has ", p_b)
            end
        end
    end
end
function test_compare(a::MJ.MarkovBiasTemperature, b::MJ.MarkovBiasTemperature)
    if a.amount != b.amount
        println(tab, "Mismatched temperature! A is ", a.amount, " while B is ", b.amount)
    end
end


DEFAULT_PRIORITY = MJ.MarkovRewritePriority_Everything()

BIG_TEST = @markovjunior 3 'R' begin
    @pragma Hi 1 3 22
    @pragma hello
    @pragma Hi "abcd"

    @rewrite R => G
    @rewrite RGB => bgw
    @rewrite 3 R => G
    @rewrite (area/2) R => _
    @rewrite (area*2) _ => R
    @rewrite RGB => [2]_[1]

    # Next op is 7
    @rewrite (1.5*area)    R[Gw]B  => b[MR]{wgb}
    @rewrite (length/2.0)  ___     => [1][3][2]  *2
    @rewrite (length*3.0)  [Rw][G] => {wgb}[M]   /1.5

    # Next op is 10
    @rewrite (4.0*length)  R=>G  %0.1       \[ +w ] #NOTE: symmetry should be ignored by parser because the rule is a single pixel
    @rewrite (2:10)        R=>G  %0.2  *3.5
    @rewrite               R=>G  %0.3  /4.1

    # Next op is 13
    @rewrite ((area*4.2):(0.5*length)) RM=>GT  \[ x ]
    @rewrite (8:(area/4.2))            RM=>GT  \[ +(1) ]
    @rewrite                           RgR=>GMG                   temperature(0.2)
    @rewrite                           RgR=>G[1]G  \[ x, -(2), 4... ] temperature(0.1)

    # Next op is 17
    @rewrite begin
        PRIORITIZE(rare)
        R => G
        R_[Bb]w => [2]_[bB]{wbR}  %(0.4:0.6)  *0.2   \[ -z..., +(8)...]
    end begin
        temperature(40.9)
    end

    # Next op is 18
    @fill 'R' uv(min=0, size=0.2)
    @fill 'b' uv(size=1, center=0) %0.2
    @fill 'w' uv(size=(0.1, 0.5), max=1) +R
    @fill 'M' pixel(min=1, max=5) -wgb %(0.1:0.9)
    @fill 'w' +gbw

    # Next op is 23
    @sequence @rewrite(10, R=>b)
    @sequence (area/2) @rewrite(10, R=>R) begin
        temperature(11.2)
    end
    @sequence begin
        @rewrite R => G
        @rewrite begin
            PRIORITIZE(earliest)
            G => B
            G => Y *2
        end temperature(0.4)
    end temperature(0.9)

    # Next op is 26
    @rewrite [ [Rw] _ B ] => [ R G w ]
    @rewrite [ R G B
               w g b
             ] => [
                _ _ _
                _ _ _
             ]
    @rewrite [
        R G B
        w g b
      ] => [
        _ _ _
        _ _ _
      ] %0.1 /2 \[
        {x, y, z}
    ]
    @rewrite(
        [ R G B
          w g b ;;;
          M T O
          [MTO] R R
        ] => [
          [1, 2, 1] [1, 2, 2] [2, 1, 1]
          [2, 1, 1] [3, 2, 2] [3, 2, 2] ;;;
          {RGB} {MTO} {wgb}
          [RGB]      G     B
        ] \[
            x[ +y, +z... ],
            (2, z)[ (+(1), -(3)), (-(3), +(1)) ]
        ]
    )

    # Next op is 30
    @rewrite w=>b field(RG)
    @rewrite w=>b field(-RG)
    @rewrite w=>b field(RG->B)
    @rewrite w=>b field(RG<-B)
    @rewrite w=>b field(-RG->B)
    # Next op is 35
    @rewrite w=>b field(RG->BR & gw)
    @rewrite w=>b field(RG<-BR & wb)
    @rewrite w=>b field(-RG->BR & bgw)
    # Next op is 38
    @rewrite w=>b begin
        temperature(5.0)
        field(R -> GB & B,
            live,
            soft=2.01,
            combo=diff
        )
        field(w -> b,
            forbidden,
            randomness=0.3,
            scale=0.5
        )
    end
end

CELL_CODE = MJ.CELL_CODE_BY_CHAR
WILDCARD = MJ.RewriteRuleCell_Wildcard()
BIG_TEST_ANSWER = MJ.MarkovAlgorithm(
    CELL_CODE['R'],
    3, 3,

    [
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            nothing,
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['b']),
                        (CELL_CODE['G'], CELL_CODE['g']),
                        (CELL_CODE['B'], CELL_CODE['w'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 1
                )
            ),
            nothing,
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            3,
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], WILDCARD)
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdByArea(0.5f0),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (WILDCARD, CELL_CODE['R'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdByArea(2.0f0),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], MJ.RewriteRuleCell_Lookup{Int}(2)),
                        (CELL_CODE['G'], WILDCARD),
                        (CELL_CODE['B'], MJ.RewriteRuleCell_Lookup{Int}(1))
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 1
                )
            ),
            nothing,
            ()
        ),

        # Next op is 7
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['b']),
                        (
                            MJ.RewriteRuleCell_Set('w', 'G'),
                            MJ.RewriteRuleCell_List((
                                CELL_CODE['R'],
                                CELL_CODE['M']
                            ))
                        ),
                        (
                            CELL_CODE['B'],
                            MJ.RewriteRuleCell_Set('w', 'g', 'b')
                        )
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 1
                )
            ),
            MJ.ThresholdByArea(1.5f0),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (WILDCARD, MJ.RewriteRuleCell_Lookup{Int}(1)),
                        (WILDCARD, MJ.RewriteRuleCell_Lookup{Int}(3)),
                        (WILDCARD, MJ.RewriteRuleCell_Lookup{Int}(2))
                    ),
                    nothing, 2.0f0,
                    MJ.GridDir[ ], 1
                )
            ),
            MJ.ThresholdByLength(0.5f0),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (
                            MJ.RewriteRuleCell_Set('w', 'R'),
                            MJ.RewriteRuleCell_Set('w', 'g', 'b')
                        ),
                        (
                            MJ.RewriteRuleCell_Set('G'),
                            MJ.RewriteRuleCell_List(tuple(
                                CELL_CODE['M']
                            ))
                        )
                    ),
                    nothing, convert(Float32, 1 / 1.5),
                    MJ.GridDir[ ], 1
                )
            ),
            MJ.ThresholdByLength(3.0f0),
            ()
        ),

        # Next op is 10
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    0.1f0, 1.0f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdByLength(4.0f0),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    0.2f0, 3.5f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdRange(2, 10),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    0.3f0, convert(Float32, 1 / 4.1),
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                )
            ),
            nothing,
            ()
        ),

        # Next op is 13
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G']),
                        (CELL_CODE['M'], CELL_CODE['T'])
                    ),
                    nothing, 1.0f0,
                    [ MJ.GridDir(1, -1), MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdRange(
                MJ.ThresholdByArea(4.2f0),
                MJ.ThresholdByLength(0.5f0)
            ),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G']),
                        (CELL_CODE['M'], CELL_CODE['T'])
                    ),
                    nothing, 1.0f0,
                    [ MJ.GridDir(1, 1) ], nothing
                )
            ),
            MJ.ThresholdRange(
                8,
                MJ.ThresholdByArea(convert(Float32, 1/4.2))
            ),
            ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G']),
                        (CELL_CODE['g'], CELL_CODE['M']),
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], (nothing, 1)
                )
            ),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(0.2)
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G']),
                        (CELL_CODE['g'], MJ.RewriteRuleCell_Lookup(1)),
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    nothing, 1.0f0,
                    [ MJ.GridDir(1, -1), MJ.GridDir(1, 1), MJ.GridDir(2, -1) ], 4
                )
            ),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(0.1)
            )
        ),

        # Next op is 17
        MJ.MarkovOpRewrite(
            MJ.MarkovRewritePriority_Rare(),
            tuple(
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], CELL_CODE['G'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                ),
                MJ.RewriteRule_Strip(
                    tuple(
                        (CELL_CODE['R'], MJ.RewriteRuleCell_Lookup{Int}(2)),
                        (WILDCARD, WILDCARD),
                        (
                            MJ.RewriteRuleCell_Set('b', 'B'),
                            MJ.RewriteRuleCell_List((
                                CELL_CODE['B'],
                                CELL_CODE['b']
                            ))
                        ),
                        (
                            CELL_CODE['w'],
                            MJ.RewriteRuleCell_Set('w', 'b', 'R')
                        )
                    ),
                    (0.4f0, 0.6f0), 0.2f0,
                    MJ.GridDir[ ], (3, 8)
                )
            ),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(40.9f0)
            )
        ),

        # Next op is 18
        MJ.MarkovOpDrawBox(
            CELL_CODE['R'],
            MJ.DrawBoxSpace.uv,
            true,
            Box1Df(min=Vec(0), size=Vec(0.2)),
            nothing, nothing
        ),
        MJ.MarkovOpDrawBox(
            CELL_CODE['b'],
            MJ.DrawBoxSpace.uv,
            true,
            Box1Df(size=Vec(1), center=Vec(0)),
            nothing, 0.2f0
        ),
        MJ.MarkovOpDrawBox(
            CELL_CODE['w'],
            MJ.DrawBoxSpace.uv,
            false,
            Box2Df(size=Vec(0.1f0, 0.5f0), max=Vec(1.0f0, 1.0f0)),
            (Val(:whitelist), MJ.CellTypeSet("R")),
            nothing
        ),
        MJ.MarkovOpDrawBox(
            CELL_CODE['M'],
            MJ.DrawBoxSpace.pixel,
            true,
            Box1Di(min=Vec(1), max=Vec(5)),
            (Val(:blacklist), MJ.CellTypeSet("wgb")),
            (0.1f0, 0.9f0)
        ),
        MJ.MarkovOpDrawBox(
            CELL_CODE['w'],
            MJ.DrawBoxSpace.uv,
            true,
            Box1Df(min=Vec(0), max=Vec(1)),
            (Val(:whitelist), MJ.CellTypeSet("gbw")),
            nothing
        ),

        # Next op is 23
        MJ.MarkovOpSequence(
            MJ.AbstractMarkovOp[
                MJ.MarkovOpRewrite(
                    DEFAULT_PRIORITY,
                    tuple(
                        MJ.RewriteRule_Strip(
                            tuple(
                                (CELL_CODE['R'], CELL_CODE['b'])
                            ),
                            nothing, 1.0f0,
                            MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                        )
                    ),
                    10,
                    ()
                )
            ],
            nothing,
            ()
        ),
        MJ.MarkovOpSequence(
            MJ.AbstractMarkovOp[
                MJ.MarkovOpRewrite(
                    DEFAULT_PRIORITY,
                    tuple(
                        MJ.RewriteRule_Strip(
                            tuple(
                                (CELL_CODE['R'], CELL_CODE['R'])
                            ),
                            nothing, 1.0f0,
                            MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                        )
                    ),
                    10,
                    ()
                )
            ],
            MJ.ThresholdByArea(0.5f0),
            tuple(
                MJ.MarkovBiasTemperature(convert(Float32, 11.2))
            )
        ),
        MJ.MarkovOpSequence(
            MJ.AbstractMarkovOp[
                MJ.MarkovOpRewrite(
                    DEFAULT_PRIORITY,
                    tuple(
                        MJ.RewriteRule_Strip(
                            tuple(
                                (CELL_CODE['R'], CELL_CODE['G'])
                            ),
                            nothing, 1.0f0,
                            MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                        )
                    ),
                    nothing,
                    ()
                ),
                MJ.MarkovOpRewrite(
                    MJ.MarkovRewritePriority_Earliest(),
                    tuple(
                        MJ.RewriteRule_Strip(
                            tuple(
                                (CELL_CODE['G'], CELL_CODE['B'])
                            ),
                            nothing, 1.0f0,
                            MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                        ),
                        MJ.RewriteRule_Strip(
                            tuple(
                                (CELL_CODE['G'], CELL_CODE['Y'])
                            ),
                            nothing, 2.0f0,
                            MJ.GridDir[ MJ.GridDir(1, 1) ], nothing
                        )
                    ),
                    nothing,
                    tuple(
                        MJ.MarkovBiasTemperature(convert(Float32, 0.4))
                    )
                )
            ],
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(0.9f0)
            )
        ),

        # Next op is 26
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_MD(
                permutedims(MJ.RewriteCell_MD{2}[ (MJ.RewriteRuleCell_Set('w', 'R'), CELL_CODE['R'])  (WILDCARD, CELL_CODE['G'])   (CELL_CODE['B'], CELL_CODE['w'])   ]),
                nothing, 1.0f0, MJ.RewriteRule_MD_Symmetry_Definition()
            )),
            nothing, ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_MD(
                permutedims(MJ.RewriteCell_MD{2}[
                    (CELL_CODE['R'], WILDCARD)    (CELL_CODE['G'], WILDCARD)    (CELL_CODE['B'], WILDCARD)
                    (CELL_CODE['w'], WILDCARD)    (CELL_CODE['g'], WILDCARD)    (CELL_CODE['b'], WILDCARD)
                ]),
                nothing, 1.0f0, MJ.RewriteRule_MD_Symmetry_Definition()
            )),
            nothing, ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_MD(
                permutedims(MJ.RewriteCell_MD{2}[
                    (CELL_CODE['R'], WILDCARD)    (CELL_CODE['G'], WILDCARD)    (CELL_CODE['B'], WILDCARD)
                    (CELL_CODE['w'], WILDCARD)    (CELL_CODE['g'], WILDCARD)    (CELL_CODE['b'], WILDCARD)
                ]),
                0.1f0,
                0.5f0,
                MJ.RewriteRule_MD_Symmetry_Definition(
                    Vector{Pair{Matrix{Int}, MJ.RewriteRule_TailSymmetry}}(),
                    Set{Int}[
                        Set([ 1, 2, 3 ])
                    ]
                )
            )),
            nothing, ()
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_MD(
                permutedims(MJ.RewriteCell_MD{3}[
                    (CELL_CODE['R'], MJ.RewriteRuleCell_Lookup((1, 2, 1)))       (CELL_CODE['G'], MJ.RewriteRuleCell_Lookup((1, 2, 2)))       (CELL_CODE['B'], MJ.RewriteRuleCell_Lookup((2, 1, 1)))
                    (CELL_CODE['w'], MJ.RewriteRuleCell_Lookup((2, 1, 1)))       (CELL_CODE['g'], MJ.RewriteRuleCell_Lookup((3, 2, 2)))       (CELL_CODE['b'], MJ.RewriteRuleCell_Lookup((3, 2, 2)))   ;;;
                    (CELL_CODE['M'], MJ.RewriteRuleCell_Set('R', 'G', 'B'))      (CELL_CODE['T'], MJ.RewriteRuleCell_Set('M', 'T', 'O'))      (CELL_CODE['O'], MJ.RewriteRuleCell_Set('w', 'g', 'b'))
                    (MJ.RewriteRuleCell_Set('M', 'T', 'O'), MJ.RewriteRuleCell_List((c->CELL_CODE[c]).(('R', 'G', 'B'))))    (CELL_CODE['R'], CELL_CODE['G'])    (CELL_CODE['R'], CELL_CODE['B'])
                ], (2, 1, 3)),
                nothing, 1.0f0,
                MJ.RewriteRule_MD_Symmetry_Definition(
                    Pair{Matrix{Int}, MJ.RewriteRule_TailSymmetry}[
                        [
                            1    2
                        ] => (nothing, 3),
                        [
                            2    1   -3
                            3    -3   1
                        ] => nothing
                    ],
                    Vector{Set{Int}}()
                )
            )),
            nothing, ()
        ),

        # Next op is 30
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, false,

                    MJ.CellTypeSet(), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet(),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, true,

                    MJ.CellTypeSet(), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet(),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, false,

                    MJ.CellTypeSet("B"), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet(),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, true,

                    MJ.CellTypeSet("B"), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet(),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, true,

                    MJ.CellTypeSet("B"), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet(),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        # Next op is 35
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, false,

                    MJ.CellTypeSet("BR"), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet("gw"),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, true,

                    MJ.CellTypeSet("BR"), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet("wb"),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("RG"),
                    MJ.BiasFieldComboMode.average,

                    false, true,

                    MJ.CellTypeSet("BR"), MJ.BiasFieldOutsideCellMode.normal,

                    MJ.CellTypeSet("bgw"),

                    0.0f0, 1.0f0, 2.0f0
                )
            )
        ),
        # Next op is 38
        MJ.MarkovOpRewrite(
            DEFAULT_PRIORITY,
            tuple(MJ.RewriteRule_Strip(tuple( (CELL_CODE['w'], CELL_CODE['b']) ),
                                       nothing, 1.0f0,
                                       MJ.GridDir[ MJ.GridDir(1, 1) ], nothing)),
            nothing,
            tuple(
                MJ.MarkovBiasTemperature(5.0),
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("R"),
                    MJ.BiasFieldComboMode.diff,

                    true, false,

                    MJ.CellTypeSet("GB"), MJ.BiasFieldOutsideCellMode.soft,

                    MJ.CellTypeSet("B"),

                    0.0f0, 1.0f0, convert(Float32, 2.01)
                ),
                MJ.MarkovBiasField(
                    MJ.CellTypeSet("w"),
                    MJ.BiasFieldComboMode.average,

                    false, true,

                    MJ.CellTypeSet("b"), MJ.BiasFieldOutsideCellMode.flippable,

                    MJ.CellTypeSet(),

                    0.3f0, 0.5f0, 2.0f0
                )
            )
        )
    ],

    [ :Hi => 1, :hello => 1, :Hi => 2 ],
    Dict(
        :Hi => Any[
            Any[ 1, 3, 22 ],
            Any[ "abcd" ]
        ],
        :hello => Any[
            Any[ ]
        ]
    ),
    Dict{Symbol, Any}()
)

@bp_check(BIG_TEST == BIG_TEST_ANSWER,
          test_compare(BIG_TEST, BIG_TEST_ANSWER),
          "INVALID result from `@markovjunior`! ",
            "Detailed printout is above this line -- A is the actual, B is the expected")

# Test dsl_string() by executing it, parsing the result, and comparing them again.
BIG_TEST_2 = markov_algo_parse(MJ.dsl_string(BIG_TEST))
@bp_check(BIG_TEST_2 == BIG_TEST_ANSWER,
          test_compare(BIG_TEST_2, BIG_TEST_ANSWER),
          "INCORRECT result from `dsl_string()`! ",
            "Detailed printout is above this line -- A is parsed from `dsl_string()`, B is the expected")

# A == B and B == C; do a sanity-check by comparing A == C.
@bp_check(BIG_TEST == BIG_TEST_2,
          test_compare(BIG_TEST, BIG_TEST_2),
          "SANITY FAIL! Transitive equality seems to be broken (a==b and b==c, but a!=b). ",
            "Detailed comparison printout is above this line -- A is from `@markovjunior` and B is parsed from `dsl_string()`")

end)()