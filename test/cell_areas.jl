# Tests for the kinds of `CellArea`, and that the built-in biases handle all of them.

(() -> begin

VI2 = MJ.CellIdx{2}
c(ch::Char) = MJ.CELL_CODE_BY_CHAR[ch]

"Collects `(local_idx, cell)` for every cell in the area"
visit_all(area) = let output = Any[ ]
    MJ.for_each_cell(area) do local_idx, cell
        push!(output, (local_idx, cell))
    end
    output
end


# CellList: visits its cells in order, with 1D local indices.
list = MJ.CellList(VI2[ VI2(3, 4), VI2(1, 2), VI2(5, 1) ])
@bp_check(visit_all(list) == [ (Int32(1), VI2(3, 4)), (Int32(2), VI2(1, 2)), (Int32(3), VI2(5, 1)) ],
          visit_all(list))
@bp_check(MJ.first_cell(list) == VI2(3, 4))

# CellMask: visits only its `true` cells, with local indices into the mask grid.
mask = MJ.CellMask(Bool[ 0 1 ; 1 1 ; 0 1 ], VI2(10, 20)) # 3x2, where Julia's first axis is X here
@bp_check(visit_all(mask) == [ (VI2(2, 1), VI2(12, 21)),
                               (VI2(1, 2), VI2(11, 22)), (VI2(2, 2), VI2(12, 22)), (VI2(3, 2), VI2(13, 22)) ],
          visit_all(mask))
@bp_check(MJ.first_cell(mask) == VI2(12, 21), MJ.first_cell(mask))

# Breaking early.
n_visited = Ref(0)
@bp_check(MJ.for_each_cell(mask, true) do local_idx, cell
    n_visited[] += 1
    return n_visited[] == 2
end)
@bp_check(n_visited[] == 2, n_visited[])
@bp_check(!MJ.for_each_cell((l, c) -> false, list, true))


###########
#  Biases

"Makes an algorithm state for the given grid, for evaluating biases directly"
make_state(grid) = MJ.AlgoState(grid, MJ.PRNG(0x1234), MJ.MarkovAllocatorHeap(),
                                Dict{Symbol, Any}(), MJ.AlgoCommsChannel(0), MarkovTickSettings(),
                                Any[ ], 0, 0, 0, 0)
"Gets the bias from a `@rewrite` written with exactly one"
parse_bias(src::String) = only(markov_algo_parse("@markovjunior 2 begin\n @rewrite b => w $src\nend").sequence[1].biases)
algo = @markovjunior 2 begin
    @rewrite b => w
end
function bias_at(bias, algo_state, area)
    state = MJ.markov_bias_initialize(bias, algo, algo_state)
    try
        return MJ.markov_bias_calculate(bias, state, algo, algo_state, area)
    finally
        MJ.markov_bias_cleanup(bias, state, algo, algo_state)
    end
end

grid = fill(c('b'), 8, 8)
grid[2, 3] = c('Y')
grid[7, 6] = c('R')
algo_state = make_state(grid)

# Every kind of area should get the same bias for the same cells.
box = MJ.BoxI{2}(min=VI2(3, 2), size=VI2(4, 3))
box_cells = vec(collect(VI2, min_inclusive(box):max_inclusive(box)))
as_region = MJ.CellRegion(box)
as_list = MJ.CellList(box_cells)
as_mask = MJ.CellMask(fill(true, 4, 3), VI2(2, 1))
as_padded_mask = MJ.CellMask(let m = fill(false, 6, 5); m[2:5, 2:4] .= true; m end, VI2(1, 0))
for combo in ("average", "min", "max", "first", "deviation", "diff")
    for src in ("field(Y, combo=$combo)", "field(-Y, combo=$combo)", "field(Y->b, combo=$combo)")
        bias = parse_bias(src)
        expected = bias_at(bias, algo_state, as_region)
        for area in (as_list, as_mask, as_padded_mask)
            actual = bias_at(bias, algo_state, area)
            @bp_check(actual == expected, src, " on ", typeof(area), ": ", actual, " vs ", expected)
        end
    end
end

# An L-shaped mask matches the equivalent list.
l_mask = MJ.CellMask(Bool[ 1 0 0 ; 1 0 0 ; 1 1 1 ], VI2(4, 4))
l_list = MJ.CellList(VI2[ cell for (_, cell) in visit_all(l_mask) ])
for src in ("field(Y)", "field(R, combo=max)", "field(Y, combo=first)")
    bias = parse_bias(src)
    @bp_check(bias_at(bias, algo_state, l_mask) == bias_at(bias, algo_state, l_list), src)
end

# A pathed field forbids areas entirely off the path, whatever their shape.
grid[1:8, 8] .= c('g')
algo_state = make_state(grid)
bias = parse_bias("field(Y->b)")
@bp_check(isnothing(bias_at(bias, algo_state, MJ.CellList(VI2[ VI2(1, 8), VI2(8, 8) ]))))
@bp_check(isnothing(bias_at(bias, algo_state, MJ.CellMask(fill(true, 8, 1), VI2(0, 7)))))
@bp_check(exists(bias_at(bias, algo_state, MJ.CellMask(fill(true, 1, 2), VI2(0, 6)))))

# Temperature doesn't care about the area.
bias = parse_bias("temperature(3)")
for area in (as_region, as_list, as_mask, l_mask)
    @bp_check(0 <= bias_at(bias, algo_state, area) <= 3)
end

end)()
