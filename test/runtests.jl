# Make sure the test is always running in the same directory and within the same project.
using Pkg
const MAIN_PROJECT_DIR = joinpath(@__DIR__, "..")
cd(MAIN_PROJECT_DIR)
Pkg.activate(".")

using MarkovJunior; const MJ = MarkovJunior
MJ.markovjunior_asserts_enabled() = true

using Bplus; @using_bplus


const DEFAULT_PRIORITY = MJ.MarkovRewritePriority_Everything()

const BIG_TEST = @markovjunior 3 'R' begin
    @pragma Hi 1 3 22
    @pragma hello

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
    @rewrite                           RM=>GT                   temperature(0.2)
    @rewrite                           RM=>GT  \[ x, -(2), 4... ] temperature(0.1)

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

    # Next op is 22
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

    # Next op is 25
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
end

const CELL_CODE = MJ.CELL_CODE_BY_CHAR
const WILDCARD = MJ.RewriteRuleCell_Wildcard()
const BIG_TEST_ANSWER = MJ.MarkovAlgorithm(
    CELL_CODE['R'],
    3, 3,

    Pair{Symbol, Vector{Any}}[
        :Hi => Any[ 1, 3, 22 ],
        :hello => Any[ ]
    ],

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
                        (CELL_CODE['M'], CELL_CODE['T'])
                    ),
                    nothing, 1.0f0,
                    MJ.GridDir[ ], 1
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
                        (CELL_CODE['M'], CELL_CODE['T'])
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

        # Next op is 13
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
            Box1Di(size=Vec(1), center=Vec(0)),
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

        # Next op is 22
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
            MJ.AbstractMarkovBias[

            ]
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
            MJ.AbstractMarkovBias[
                MJ.MarkovBiasTemperature(convert(Float32, 11.2))
            ]
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
            MJ.AbstractMarkovBias[
                MJ.MarkovBiasTemperature(0.9f0)
            ]
        ),

        # Next op is 25
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
        )
    ]
)

function test_compare(a::MJ.MarkovAlgorithm, b::MJ.MarkovAlgorithm, tab::String = "")
    println("Comparing two algorithms...")
    if length(a.sequence) != length(b.sequence)
        println(tab, "\tA has ", length(a.sequence), " elements while B has ", length(b.sequence), "!")
    else
        for (oa, ob, i) in zip(a.sequence, b.sequence, 1:length(a.sequence))
            println(tab, "\tOp ", i, ", ", typeof(oa), "(", typeof(ob), ")...")
            test_compare(oa, ob, "$tab\t\t")
        end
    end
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
            println(tab, "Bias ", i, ", ", MJ.dsl_string(ba), "   (", MJ.dsl_string(bb), ")")
            if ba != bb
                println(tab, "\tMismatch!")
            end
        end
        if a.biases != b.biases
            println(tab, "Bias list Mismatch!")
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
            println(tab, "Bias ", i, ", ", MJ.dsl_string(ba), "   (", MJ.dsl_string(bb), ")")
            if ba != bb
                println(tab, "\tMismatch!")
            end
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

@bp_check(BIG_TEST == BIG_TEST_ANSWER,
          test_compare(BIG_TEST, BIG_TEST_ANSWER),
          "INVALID result from `@markovjunior`! ",
            "Detailed printout is above this line -- A is the actual, B is the expected")
const BIG_TEST_2 = MJ.parse_markovjunior(MJ.dsl_string(BIG_TEST))
@bp_check(BIG_TEST_2 == BIG_TEST_ANSWER,
          test_compare(BIG_TEST_2, BIG_TEST_ANSWER),
          "INCORRECT result from `dsl_string()`! ",
            "Detailed printout is above this line -- A is parsed from `dsl_string()`, B is the expected")
# Do a sanity-check:
@bp_check(BIG_TEST == BIG_TEST_2,
          test_compare(BIG_TEST, BIG_TEST_2),
          "SANITY FAIL! Transitive equality seems to be broken (a==b and b==c, but a!=b). ",
            "Detailed comparison printout is above this line -- A is from `@markovjunior` and B is parsed from `dsl_string()`")


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
# FOR NEW_TESTS:
# MJ.log_md_symmetry_logic() = true


println("\n\nTests passed!\n")