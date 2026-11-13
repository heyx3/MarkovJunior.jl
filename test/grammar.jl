# Tests for `@grammar`, the split-grammar Op.
# Uses `CountingBias` from scale.jl, so it must be included after that file.

(() -> begin

c(ch::Char) = MJ.CELL_CODE_BY_CHAR[ch]
cells_of(grid, ch::Char) = Set(Tuple.(findall(==(c(ch)), grid)))
neighbors(x, y) = Set([ (x-1, y), (x+1, y), (x, y-1), (x, y+1) ])
str(grid::Vector{UInt8}) = String([ MJ.CELL_TYPES[v + 1].char for v in grid ])
"How many lines were logged"
n_logged(err::String) = count(!isempty, split(err, '\n'))

function check_roundtrip(mj::MarkovAlgorithm)
    test_parsed_mj = markov_algo_parse(markov_algo_to_string(mj))
    @bp_check(mj == test_parsed_mj,
              "Parsing failed!\n\nSource:\n", mj,
              "\nDest:\n", test_parsed_mj,
              "\nSource as string: ", markov_algo_to_string(mj))
end

"Runs the given algorithm to completion on the given grid, skipping all tagged events"
function run_on(mj::MarkovAlgorithm, initial::Array{UInt8}; seed = 0xaabbcc)
    check_roundtrip(mj)
    state_ch = markov_algo_run(mj, initial,
                               MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                               seeds=seed,
                               allocator=MJ.MarkovAllocatorHeapReused())
    return markov_algo_complete(copy, state_ch)
end
"Runs a grammar (given as the inside of a `@grammar` block) on a 1D grid, returning the grid as a string"
run_1D(grammar::String, len::Int; initial = 'b', seed = 0xaabbcc) = str(run_on(
    markov_algo_parse("@markovjunior 1 begin\n @grammar begin\n $grammar\n end\nend"),
    fill(c(initial), len);
    seed=seed
))

"Checks that parsing fails (but not the error message, which may change)"
function check_parse_fails(src::String)
    try
        markov_algo_parse(src)
    catch e
        return nothing
    end
    error("Should have failed to parse:\n", src)
end
check_grammar_parse_fails(grammar::String) = check_parse_fails(
    "@markovjunior begin\n @grammar begin\n $grammar\n end\nend"
)


#########
# split()

ABC = "\n A => 'R'\n B => 'G'\n C => 'Y'"

# Pixel sizes come first, then weights share what's left.
s = run_1D("ROOT => split(x, A(2), B(w=1), C(w=3))" * ABC, 10)
@bp_check(s == "RRGGYYYYYY", s)
# A negative axis lays the children out in reverse.
s = run_1D("ROOT => split(-x, A(2), B(w=1), C(w=3))" * ABC, 10)
@bp_check(s == "YYYYYYGGRR", s)
# Children without a size have a weight of 1.
s = run_1D("ROOT => split(x, A(2), B, C)" * ABC, 10)
@bp_check(s == "RRGGGGYYYY", s)
s = run_1D("ROOT => split(x, A(2), B(), C)" * ABC, 10)
@bp_check(s == "RRGGGGYYYY", s)
# Without any weights, pixel sizes stretch to fill the space.
s = run_1D("ROOT => split(x, A(1), B(3))" * ABC, 8)
@bp_check(s == "RRGGGGGG", s)
# `len` is the parent's length along the axis.
s = run_1D("ROOT => split(x, A(0.25len), B)" * ABC, 8)
@bp_check(s == "RRGGGGGG", s)
# A pixel size can be any calculation; `~` is bitwise-not, so this is a 2-pixel size.
s = run_1D("ROOT => split(x, A(~(-3)), B)" * ABC, 5)
@bp_check(s == "RRGGG", s)
# `_` children leave their part untouched.
s = run_1D("ROOT => split(x, A(2), _(w=1), A(2))" * ABC, 7)
@bp_check(s == "RRbbbRR", s)
# Uneven sizes are rounded so the children always fill the parent.
for seed in 1:10
    s = run_1D("ROOT => split(x, A, B, C)" * ABC, 10; seed=seed)
    @bp_check(sort([ count(==(ch), s) for ch in "RGY" ]) == [ 3, 3, 4 ], s)
    @bp_check(issorted(s; by=ch -> findfirst(ch, "RGY")), s)
end
# Random sizes.
for seed in 1:10
    s = run_1D("ROOT => split(x, A(1:3), B)" * ABC, 10; seed=seed)
    @bp_check(count(==('R'), s) in 1:3, s)
end

# Splitting a 2D grid along each axis.
grid = run_on(@markovjunior(begin
    @grammar begin
        ROOT => split(y, Bottom(1), Top(w=1))
        Bottom => 'R'
        Top => split(size_max_axis(), A(2), B)
        A => 'G'
        B => 'Y'
    end
end), fill(c('b'), 4, 7))
@bp_check(cells_of(grid, 'R') == Set((x, 1) for x in 1:4), grid)
@bp_check(cells_of(grid, 'G') == Set((x, y) for x in 1:4, y in 2:3), grid)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 1:4, y in 4:7), grid)


##############
# Alternatives

# An alternative whose pixel sizes don't fit is skipped.
s = run_1D("ROOT => split(x, A(5), B(5)) *100\n ROOT => C" * ABC, 8)
@bp_check(s == "YYYYYYYY", s)
# If nothing fits, the region is untouched.
s = run_1D("ROOT => split(x, A(5), B(5))" * ABC, 8)
@bp_check(s == "bbbbbbbb", s)
# Guards.
@bp_check(run_1D("ROOT(size(x) > 5) => A\n ROOT(size(x) <= 5) => B" * ABC, 6) == "RRRRRR")
@bp_check(run_1D("ROOT(size(x) > 5) => A\n ROOT(size(x) <= 5) => B" * ABC, 4) == "GGGG")
@bp_check(run_1D("ROOT(size_max() == grid_size(1)) => split(x, ROOT, ROOT)\n ROOT(size_max() < grid_size(1)) => A" * ABC, 4) == "RRRR")
# `depth` counts the expansions.
@bp_check(run_1D("ROOT => T\n T(depth == 1) => A\n T(depth != 1) => B" * ABC, 3) == "RRR")
# Zero weights are never picked.
@bp_check(run_1D("ROOT => A *0\n ROOT => B" * ABC, 3) == "GGG")
@bp_check(run_1D("ROOT => A *(depth)\n ROOT => B" * ABC, 3) == "GGG")
# Weights bias the choice.
for (weight, expected) in (("*3", 0.75), ("/3", 0.25))
    s = run_1D("ROOT => repeat(x, T(1))\n T => A $weight\n T => B" * ABC, 400)
    frac = count(==('R'), s) / length(s)
    @bp_check(abs(frac - expected) < 0.1, weight, ": ", frac)
end
# Recursion stops once regions get too small.
s = run_1D("ROOT(size(x) >= 2) => split(x, A(1), ROOT)\n ROOT(size(x) < 2) => B" * ABC, 5)
@bp_check(s == "RRRRG", s)
# `ROOT` covers the whole grid wherever its rules are, and other symbols can refer back to it.
@bp_check((s = run_1D(ABC * "\n ROOT => split(x, A(1), B)", 3)) == "RGG", s)
s = run_1D("ROOT(size(x) >= 2) => split(x, A(1), T)\n T => ROOT\n ROOT(size(x) < 2) => B" * ABC, 5)
@bp_check(s == "RRRRG", s)


##########
# repeat()

# As many copies as fit, stretched to fill the region.
for seed in 1:5
    s = run_1D("ROOT => repeat(x, T(3))\n T => split(x, A(1), B)" * ABC, 10; seed=seed)
    @bp_check(count(==('R'), s) == 3, s)
    @bp_check(s[1] == 'R', s)
end
# If not even one copy fits, the alternative doesn't fit.
@bp_check(run_1D("ROOT => repeat(x, A(11))\n ROOT => B" * ABC, 10) == "G"^10)
# Separators go between the copies, and keep their size while the copies stretch.
for seed in 1:5
    s = run_1D("ROOT => repeat(x, A(3), B(1))" * ABC, 12; seed=seed)
    @bp_check(occursin(r"^RRRR?GRRRR?GRRRR?$", s) && (count(==('R'), s) == 10), s)
end
# The leftover space can go to an extra child instead.
@bp_check((s = run_1D("ROOT => repeat(x, A(3), B(1); leftover=after(C))" * ABC, 12)) == "RRRGRRRGRRRY", s)
@bp_check((s = run_1D("ROOT => repeat(x, A(3), B(1); leftover=before(_))" * ABC, 13)) == "bbRRRGRRRGRRR", s)
@bp_check((s = run_1D("ROOT => repeat(x, A(3), B(1); leftover=around(C))" * ABC, 13)) == "YRRRGRRRGRRRY", s)
@bp_check((s = run_1D("ROOT => repeat(x, A(3); leftover=after(C))" * ABC, 10)) == "RRRRRRRRRY", s)
for seed in 1:5
    s = run_1D("ROOT => repeat(x, A(3), B(1); leftover=around(C))" * ABC, 14; seed=seed)
    @bp_check(occursin(r"^YY?RRRGRRRGRRRYY?$", s) && (count(==('Y'), s) == 3), s)
end
# Or to the separators.
for seed in 1:5
    s = run_1D("ROOT => repeat(x, A(1), B(2); leftover=spread)" * ABC, 15; seed=seed)
    @bp_check(occursin(r"^RGGG?RGGG?RGGG?RGGG?R$", s) && (count(==('G'), s) == 10), s)
end
# With only one copy, there's nothing to spread, so the copy stretches.
@bp_check((s = run_1D("ROOT => repeat(x, A(3), B(1); leftover=spread)" * ABC, 5)) == "RRRRR", s)
# `each()` evaluates a size again for every copy.
copy_lengths(s) = [ length(m.match) for m in eachmatch(r"RY*", s) ]
for (size_expr, expect_varied) in (("2:5", false), ("each(2:5)", true))
    varied = false
    for seed in 1:10
        s = run_1D("ROOT => repeat(x, T($size_expr), B(1))\n T => split(x, A(1), C)" * ABC, 40; seed=seed)
        lengths = copy_lengths(s)
        @bp_check(sum(lengths) + length(lengths) - 1 == 40, s)
        varied |= (maximum(lengths) - minimum(lengths)) > 1
    end
    @bp_check(varied == expect_varied, size_expr)
end


#########
# inset()

INSET = "\n W => 'N'\n C => 'Y'"
run_2D(grammar, size...) = run_on(
    markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n $grammar\n end\nend"),
    fill(c('b'), size...)
)

grid = run_2D("ROOT => inset(1, W, C)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 2:4, y in 2:5), grid)
@bp_check(count(==(c('N')), grid) == 18, grid)
# Per-axis amounts; missing axes use 0.
grid = run_2D("ROOT => inset([1], W, C)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'N') == Set((x, y) for x in (1, 5), y in 1:6), grid)
grid = run_2D("ROOT => inset([0, 2], W, C)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'N') == Set((x, y) for x in 1:5, y in (1, 2, 5, 6)), grid)
# The core can't be empty.
grid = run_2D("ROOT => inset(3, W, C)\n ROOT => _" * INSET, 5, 6)
@bp_check(all(==(c('b')), grid), grid)
# One box per face; the corners go to the faces along the lowest axis.
grid = run_2D("ROOT => inset(1, W, _)\n W => split(size_max_axis(), A(1), _)\n A => 'R'", 5, 6)
@bp_check(cells_of(grid, 'R') == Set([ (1, 1), (5, 1), (2, 1), (2, 6) ]), grid)
# Per-axis amounts can be different on the low and high sides.
grid = run_2D("ROOT => inset([[1, 2]], W, C)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'N') == Set((x, y) for x in (1, 4, 5), y in 1:6), grid)
grid = run_2D("ROOT => inset([1, [0, 2]], W, C)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 2:4, y in 1:4), grid)
# Faces can be covered by other symbols.
grid = run_2D("ROOT => inset(1, W, C, -y => F, y => _)\n F => 'R'" * INSET, 5, 6)
@bp_check(cells_of(grid, 'N') == Set((x, y) for x in (1, 5), y in 1:6), grid)
@bp_check(cells_of(grid, 'R') == Set((x, 1) for x in 2:4), grid)
@bp_check(cells_of(grid, 'b') == Set((x, 6) for x in 2:4), grid)
grid = run_2D("ROOT => inset(1, _, C, 1 => F)\n F => 'R'" * INSET, 5, 6)
@bp_check(cells_of(grid, 'R') == Set((5, y) for y in 1:6), grid)


#########
# place()

# A box within the region; the rest is split into one box per face.
grid = run_2D("ROOT => place(pixel(min=[1, 4], size=[2, 3]), C, W)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 1:2, y in 4:6), grid)
@bp_check(count(==(c('N')), grid) == 24, grid)
grid = run_2D("ROOT => place(pixel(max=[5, 2], size=2), C, W)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 4:5, y in 1:2), grid)
grid = run_2D("ROOT => place(pixel(center=3, size=1), C, W, x => F)\n F => 'R'" * INSET, 5, 5)
@bp_check(cells_of(grid, 'Y') == Set([ (3, 3) ]), grid)
@bp_check(cells_of(grid, 'R') == Set((x, y) for x in 4:5, y in 1:5), grid)
# UV boxes are relative to the region's size.
grid = run_2D("ROOT => place(uv(center=0.5, size=0.5), C, W)" * INSET, 8, 4)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 3:6, y in 2:3), grid)
# Axes past the end of a list cover the whole region, leaving no outer boxes along them.
grid = run_2D("ROOT => place(uv(min=[0.5], max=[1]), C, W)" * INSET, 6, 4)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 4:6, y in 1:4), grid)
@bp_check(cells_of(grid, 'N') == Set((x, y) for x in 1:3, y in 1:4), grid)
# Values are calculations, so they can be random or depend on the region.
# A single value is shared by every axis, unless it's wrapped in `each()`.
for (min_expr, expected) in (("1:3", Set((i, i) for i in 1:3)), ("each(1:3)", Set((x, y) for x in 1:3, y in 1:3)))
    positions = Set{Tuple{Int, Int}}()
    for seed in 1:30
        grid = run_on(markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n ROOT => place(pixel(min=$min_expr, size=2), C, W)$INSET\n end\nend"),
                      fill(c('b'), 4, 4); seed=seed)
        ys = cells_of(grid, 'Y')
        @bp_check(length(ys) == 4, grid)
        push!(positions, minimum(ys))
    end
    @bp_check(positions == expected, min_expr, ": ", positions)
end
grid = run_2D("ROOT => place(pixel(min=2, size=size_min() - 2), C, W)" * INSET, 6, 4)
@bp_check(cells_of(grid, 'Y') == Set((x, y) for x in 2:3, y in 2:3), grid)
# Boxes are clamped to the region.
grid = run_2D("ROOT => place(pixel(min=[-1], size=[3]), C, W)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set((1, y) for y in 1:6), grid)
# `texel` is one pixel in UV space; `center` is the middle of the axis in the box's own space.
grid = run_2D("ROOT => place(uv(center=0.5, size=texel), C, W)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set([ (3, 4) ]), grid)
grid = run_2D("ROOT => place(pixel(center=center, size=1), C, W)" * INSET, 5, 6)
@bp_check(cells_of(grid, 'Y') == Set([ (3, 3) ]), grid)
# If the box rounds to nothing, the alternative doesn't fit.
grid = run_2D("ROOT => place(uv(center=0.5, size=0), C, W)
 ROOT => _" * INSET, 5, 6)
@bp_check(all(==(c('b')), grid), grid)
grid = run_2D("ROOT => place(pixel(min=[9], size=[1]), C, W)
 ROOT => _" * INSET, 5, 6)
@bp_check(all(==(c('b')), grid), grid)


##########
# mirror()

# Pixel mirrors copy the source (the low half, rounded up) onto the rest.
@bp_check((s = run_1D("ROOT => mirror(x, H)\n H => split(x, A(1), B(1), C)" * ABC, 8)) == "RGYYYYGR", s)
@bp_check((s = run_1D("ROOT => mirror(x, H)\n H => split(x, A(1), B(1), C)" * ABC, 7)) == "RGYYYGR", s)
# A negative axis uses the high half as the source.
@bp_check((s = run_1D("ROOT => mirror(-x, H)\n H => split(x, A(1), B(1), C)" * ABC, 7)) == "YYGRGYY", s)
# The copy happens after the source's Ops run.
for seed in 1:5
    grid = run_on(markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n ROOT => mirror(x, H)\n H => begin\n @rewrite (2) b => R\n end\n end\nend"),
                  fill(c('b'), 8, 5); seed=seed)
    @bp_check(grid == reverse(grid, dims=1), grid)
    @bp_check(count(==(c('R')), grid) == 4, grid)
end
# Mirroring along several axes.
for seed in 1:5
    grid = run_on(markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n ROOT => mirror([x, -y], H)\n H => begin\n @rewrite (2) b => R\n end\n end\nend"),
                  fill(c('b'), 8, 6); seed=seed)
    @bp_check((grid == reverse(grid, dims=1)) && (grid == reverse(grid, dims=2)), grid)
    @bp_check(count(==(c('R')), grid) == 8, grid)
end
# A source of `_` mirrors whatever is already there.
let initial = fill(c('b'), 5, 2)
    initial[1, 1] = c('R')
    initial[2, 2] = c('G')
    grid = run_on(markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n ROOT => mirror(x, _)\n end\nend"), initial)
    @bp_check(cells_of(grid, 'R') == Set([ (1, 1), (5, 1) ]) && cells_of(grid, 'G') == Set([ (2, 2), (4, 2) ]), grid)
end

# Layout mirrors copy the boxes, which then run independently.
any_asymmetric = false
for seed in 1:10
    grid = run_on(markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n ROOT => mirror(x, H, layout)\n H => split(x, C(2), Q)\n Q => begin\n @rewrite (2) b => R\n end$INSET\n end\nend"),
                  fill(c('b'), 10, 4); seed=seed)
    @bp_check(cells_of(grid, 'Y') == Set((x, y) for x in (1, 2, 9, 10), y in 1:4), grid)
    @bp_check(count(==(c('R')), grid[1:5, :]) == 2 && count(==(c('R')), grid[6:10, :]) == 2, grid)
    any_asymmetric |= (grid != reverse(grid, dims=1))
end
@bp_check(any_asymmetric)
# Blocks get reflected.
grid = run_2D("ROOT => mirror(x, H, layout)\n H => [ R G ]", 8, 2)
@bp_check(grid == c.([ 'R' 'R' ; 'R' 'R' ; 'G' 'G' ; 'G' 'G' ; 'G' 'G' ; 'G' 'G' ; 'R' 'R' ; 'R' 'R' ]), grid)
# Boxes are cut off where they would overlap the source's middle slice.
@bp_check((s = run_1D("ROOT => mirror(x, H, layout)\n H => split(x, A(2), B)" * ABC, 5)) == "RRGRR", s)
# Pixel mirrors within layout mirrors.
any_asymmetric = false
for seed in 1:10
    grid = run_on(markov_algo_parse("@markovjunior 2 begin\n @grammar begin\n ROOT => mirror(x, H, layout)\n H => mirror(y, Q)\n Q => begin\n @rewrite (2) b => R\n end\n end\nend"),
                  fill(c('b'), 8, 6); seed=seed)
    @bp_check(grid == reverse(grid, dims=2), grid)
    any_asymmetric |= (grid != reverse(grid, dims=1))
end
@bp_check(any_asymmetric)


#####
# 3D

run_3D(grammar, size...) = run_on(
    markov_algo_parse("@markovjunior 3 begin\n @grammar begin\n $grammar\n end\nend"),
    fill(c('b'), size...)
)

grid = run_3D("ROOT => inset(1, W, C, -z => F, z => R)\n F => 'G'\n R => 'R'" * INSET, 4, 4, 4)
@bp_check(cells_of(grid, 'Y') == Set((x, y, z) for x in 2:3, y in 2:3, z in 2:3), grid)
@bp_check(cells_of(grid, 'G') == Set((x, y, 1) for x in 2:3, y in 2:3), grid)
@bp_check(cells_of(grid, 'R') == Set((x, y, 4) for x in 2:3, y in 2:3), grid)
@bp_check(count(==(c('N')), grid) == 48, grid)
grid = run_3D("ROOT => place(pixel(min=[1, 1, 3], size=2), C, W)" * INSET, 4, 4, 4)
@bp_check(cells_of(grid, 'Y') == Set((x, y, z) for x in 1:2, y in 1:2, z in 3:4), grid)
grid = run_3D("ROOT => repeat(z, C(1), W(1))" * INSET, 2, 2, 5)
@bp_check(all(grid[:, :, z] == fill(c(isodd(z) ? 'Y' : 'N'), 2, 2) for z in 1:5), grid)
grid = run_3D("ROOT => split(y, C(1), W)" * INSET, 2, 3, 2)
@bp_check(cells_of(grid, 'Y') == Set((x, 1, z) for x in 1:2, z in 1:2), grid)
grid = run_3D("ROOT => mirror(z, H)\n H => split(z, C(1), W)" * INSET, 2, 2, 6)
@bp_check(all(grid[:, :, z] == fill(c(z in (1, 6) ? 'Y' : 'N'), 2, 2) for z in 1:6), grid)


########
# Leaves

# Blocks are stretched over the region.
@bp_check(run_1D("ROOT => [ R G ]", 4) == "RRGG")
@bp_check(run_1D("ROOT => [ R G ]", 5) == "RRRGG")
@bp_check(run_1D("ROOT => [ _ R ]", 4; initial='Y') == "YYRR")
for seed in 1:5
    s = run_1D("ROOT => [ {RG} {RG} ]", 20; seed=seed)
    @bp_check(all(in("RG"), s) && ('R' in s) && ('G' in s), s)
end
grid = run_2D("ROOT => [ G ; N ]", 3, 4)
@bp_check(cells_of(grid, 'G') == Set((x, y) for x in 1:3, y in 1:2), grid)
@bp_check(cells_of(grid, 'N') == Set((x, y) for x in 1:3, y in 3:4), grid)

# Op leaves run within their region.
grid = run_on(@markovjunior(begin
    @grammar begin
        ROOT => split(x, A(3), B)
        A => @rewrite b => w
        B => begin
            @fill 'R' uv(min=0, max=1)
            @rewrite 1 R => Y
        end
    end
end), fill(c('b'), 8, 2))
@bp_check(cells_of(grid, 'w') == Set((x, y) for x in 1:3, y in 1:2), grid)
@bp_check(count(==(c('Y')), grid) == 1, grid)
@bp_check(count(==(c('R')), grid) == 9, grid)

# A grammar inside a scope only covers the scope.
grid = run_on(@markovjunior(begin
    @within_box pixel(min=3, size=4) begin
        @grammar begin
            ROOT => split(x, A(1), B)
            A => 'R'
            B => 'G'
        end
    end
end), fill(c('b'), 8, 8))
@bp_check(cells_of(grid, 'R') == Set((3, y) for y in 3:6), grid)
@bp_check(cells_of(grid, 'G') == Set((x, y) for x in 4:6, y in 3:6), grid)

# The example from the docs.
grid = run_on(@markovjunior(begin
    @grammar begin
        ROOT => split(y, Ground(3), Floors(w=1), Roof(0.15len))
        Floors => repeat(y, Floor(6))
        Floor => inset(1, Wall, Interior)
        Wall => split(size_max_axis(), Solid(w=1), Door(1), Solid(w=1))  /2
        Wall => Solid
        Solid => 'N'
        Door => 'b'
        Interior(size_max() > 8) => split(size_max_axis(), Room(w=(1:3)), Interior(w=2))
        Interior => Room
        Room => 'Y'
        Room => begin
            @fill 'R' uv(center=0.5, size=texel)
        end  *3
        Ground => [ G ; N ]
        Roof => _
    end
end), fill(c('b'), 30, 40))
@bp_check(all(==(c('G')), grid[:, 1:2]), grid)
@bp_check(all(==(c('N')), grid[:, 3]), grid)
@bp_check(all(==(c('b')), grid[:, 35:40]), grid)
@bp_check(count(==(c('N')), grid[1, 4:34]) >= 25, grid)


########
# Biases

# Op leaves can inherit biases, getting separate states.
COUNTING_BIAS_N_INITIALIZED[] = 0
grid = run_on(@markovjunior(begin
    @sequence begin
        @grammar inherit_biases begin
            ROOT => split(x, A, A)
            A => @rewrite 1 b => R
        end
    end counting()
end), fill(c('b'), 8, 8))
@bp_check(count(==(c('R')), grid) == 2, grid)
@bp_check(COUNTING_BIAS_N_INITIALIZED[] == 3, COUNTING_BIAS_N_INITIALIZED[])
@bp_check(isempty(COUNTING_BIAS_LIVE_STATES), COUNTING_BIAS_LIVE_STATES)

grid_with_Y = let g = fill(c('b'), 8, 8)
    g[7, 7] = c('Y')
    g
end
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @sequence begin
            @grammar inherit_biases begin
                ROOT => split(x, _(4), A)
                A => @rewrite 1 b => R
            end
        end field(Y)
    end), grid_with_Y; seed=seed)
    @bp_check(only(cells_of(grid, 'R')) in neighbors(7, 7), seed, ": ", grid)
end

# Live inherited biases hear about painted leaves.
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @sequence begin
            @grammar begin
                ROOT => split(x, Left(1), _, Right(1))
                Left => split(y, Corner(1), _)
                Corner => 'Y'
                Right => split(y, _, Covered(1))
                Covered => 'w'
            end
            @rewrite 1 b => M
        end field(Y, live)
    end), let g = fill(c('b'), 8, 8); g[8, 8] = c('Y'); g end; seed=seed)
    @bp_check(only(cells_of(grid, 'M')) in neighbors(1, 1), seed, ": ", grid)
end


################
# Runtime errors

# Rules that recurse forever stop at the maximum depth, painting those regions magenta and logging once.
for (rules, expected) in (("ROOT => L1\n L1 => ROOT", "MMMM"),                          # Rename cycle
                          ("ROOT => split(x, ROOT)", "MMMM"),                           # Never shrinks
                          ("ROOT => split(x, A(2), L1)\n L1 => L2\n L2 => L1", "RRMM"),  # Only one branch loops
                          ("ROOT => repeat(x, ROOT(1))", "MMMM"))
    result = Ref{String}()
    err = @capture_err (result[] = run_1D(rules * ABC, 4))
    @bp_check(result[] == expected, rules, " gave ", result[])
    @bp_check(n_logged(err) == 1, "Expected one logged error: ", err)
end
# Deep recursion within the limit is fine.
let s = run_1D("ROOT(size(x) > 1) => split(x, A(1), ROOT(w=1))\n ROOT(size(x) <= 1) => B" * ABC, 300)
    @bp_check(s == ("R"^299) * "G", s)
end

# Errors are logged, and a fallback value is used.
result = Ref{String}()
err = @capture_err (result[] = run_1D("ROOT(size(x)) => A\n ROOT => B" * ABC, 4))
@bp_check(result[] == "GGGG", result[])
@bp_check(n_logged(err) > 0, "Nothing was logged")
err = @capture_err (result[] = run_1D("ROOT => split(ndims + 1, A)" * ABC, 4))
@bp_check(result[] == "bbbb", result[])
@bp_check(n_logged(err) > 0, "Nothing was logged")
err = @capture_err (result[] = run_1D("ROOT => split(x, A(size(3)), B)" * ABC, 4))
@bp_check(result[] == "GGGG", result[])
@bp_check(n_logged(err) > 0, "Nothing was logged")


##############
# Parsing

# Constant axes and blocks put a floor on the dimensionality.
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => split(z, A)\n A => 'R'\n end\nend").min_dimension == 3)
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => split(-y, A)\n A => 'R'\n end\nend").min_dimension == 2)
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => [ R ; G ]\n end\nend").min_dimension == 2)
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => split(size_max_axis(), A)\n A => 'R'\n end\nend").min_dimension == 1)
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => inset(1, A, A, -z => A)\n A => 'R'\n end\nend").min_dimension == 3)
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => mirror([x, -z], A)\n A => 'R'\n end\nend").min_dimension == 3)
@bp_check(markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => place(pixel(min=[1, 1, 1], size=1), A, A)\n A => 'R'\n end\nend").min_dimension == 3)

check_roundtrip(markov_algo_parse("""@markovjunior begin
    @grammar inherit_biases begin
        ROOT((size_min() > 3) && (depth < 4)) => split(-(1:2), A(curve(1:3, 2.0)), _(w=(1 + depth)), ROOT) *(2 * depth)
        ROOT => inset([1, 2.5], A, _) /3
        A => repeat(y, ROOT(len / 3))
        A => repeat(-x, ROOT(each(1:3)), _(1); leftover=around(A))
        A => repeat(x, ROOT(2), A(each(0:1)); leftover=spread)
        A => inset([[1, 2], 0, [0:1, 1]], A, ROOT, -y => _, 3 => A)
        A => place(uv(center=[0.5, 0.2:0.8, 1], size=0.3), A, ROOT, x => _)
        A => place(pixel(min=1, max=size_min() - 1), _, ROOT)
        A => mirror([x, -y], ROOT, layout)
        A => mirror(-z, _)
        A => [ R ; {GB} ;; _ ; Y ]
    end
end"""))

# Lists can be written with brackets or parentheses.
for (brackets, parens) in (("inset([1, [0, 2]], A, A)", "inset((1, (0, 2)), A, A)"),
                           ("place(uv(min=[0, 0.5], size=0.5), A, A)", "place(uv(min=(0, 0.5), size=0.5), A, A)"),
                           ("mirror([x, -y], A)", "mirror((x, -y), A)"))
    parse_grammar(production) = markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => $production\n A => 'R'\n end\nend")
    @bp_check(parse_grammar(brackets) == parse_grammar(parens), parens)
end

check_grammar_parse_fails("ROOT => split(x, A, B)\n A => 'R'")
check_grammar_parse_fails("_ => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A(w=1))\n A => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A())\n A => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A(1), B)\n A => 'R'\n B => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A(1); leftover=spread)\n A => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A(1); leftover=sideways)\n A => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A(1); fill=1)\n A => 'R'")
check_grammar_parse_fails("ROOT => repeat(x, A(1), leftover=stretch, leftover=stretch)\n A => 'R'")
check_grammar_parse_fails("ROOT => split(x, A(each(1)))\n A => 'R'")
check_grammar_parse_fails("ROOT => inset(1, A, A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => inset([[1, 2, 3]], A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => inset(1, A, A, 0 => A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(uv(min=0, size=1), A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(A, A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(foo(min=0, size=1), A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(uv(min=0), A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(uv(min=0, scale=1), A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(uv(min=[0, 1], size=[1]), A, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => place(uv(min=0, size=1), A, A; at=1)\n A => 'R'")
check_grammar_parse_fails("ROOT => mirror(x, A, pixels)\n A => 'R'")
check_grammar_parse_fails("ROOT => mirror([], A)\n A => 'R'")
check_grammar_parse_fails("ROOT => split(0, A)\n A => 'R'")
check_grammar_parse_fails("ROOT => split(1.5, A)\n A => 'R'")
check_grammar_parse_fails("ROOT(len > 3) => 'R'")
check_grammar_parse_fails("ROOT => split(x, A(nope))\n A => 'R'")
check_grammar_parse_fails("ROOT => @upscale (2)")
check_grammar_parse_fails("ROOT => 3")
check_grammar_parse_fails("ROOT => '%'")
check_grammar_parse_fails("ROOT + 1")
check_grammar_parse_fails("ROOT(a, b) => 'R'")
check_grammar_parse_fails("A => 'R'")
check_parse_fails("@markovjunior begin\n @grammar begin\n end\nend")
check_parse_fails("@markovjunior begin\n @grammar inherit_biases\nend")


#####################
# User functions (@calc)

CALC_SRC = """@markovjunior 1 begin
    @calc half(a) = a / 2
    @calc stripe() = 2
    @calc long() = grid_size_max_axis()
    @grammar begin
        ROOT => split(long(), A(half(len)), B(stripe()), C)
        A => 'R'
        B => 'G'
        C => 'Y'
    end
end"""
calc_mj = markov_algo_parse(CALC_SRC)
@bp_check(collect(keys(calc_mj.calc_macros)) == [ :half, :stripe, :long ]) # Kept in order
@bp_check(calc_mj.calc_macros[:half] == MJ.CalcInterpreter.CalcMacro([ :a ], :(a / 2)))
run_calc(mj) = str(run_on(mj, fill(c('b'), 10)))
@bp_check(run_calc(calc_mj) == "RRRRRGGYYY", run_calc(calc_mj)) # Also checks the round-trip

# Redefining a function changes the next runs, without re-parsing.
push!(calc_mj.calc_macros, calc_macro_parse("stripe() = 4"))
@bp_check(run_calc(calc_mj) == "RRRRRGGGGY", run_calc(calc_mj))
push!(calc_mj.calc_macros, calc_macro_parse(:( half(a) = a / 5 )))
@bp_check(run_calc(calc_mj) == "RRGGGGYYYY", run_calc(calc_mj))
# Functions can be added that the source never defined, and used by functions that it did.
push!(calc_mj.calc_macros, calc_macro_parse("stripe() = extra() + 1"))
push!(calc_mj.calc_macros, calc_macro_parse("extra() = depth"))
@bp_check(run_calc(calc_mj) == "RRGYYYYYYY", run_calc(calc_mj))

# A missing or broken function is a runtime error, logged like any other.
delete!(calc_mj.calc_macros, :stripe)
let result = Ref{String}()
    err = @capture_err (result[] = run_calc(calc_mj))
    @bp_check(result[] == "RRYYYYYYYY", result[])
    @bp_check(n_logged(err) > 0, "Nothing was logged")
end
push!(calc_mj.calc_macros, calc_macro_parse("stripe() = size(9)"))
let result = Ref{String}()
    err = @capture_err (result[] = run_calc(calc_mj))
    @bp_check(result[] == "RRYYYYYYYY", result[])
    @bp_check(n_logged(err) > 0, "Nothing was logged")
end
# The function's body sees the environment of the calculation calling it, e.g. `len` only in child sizes.
let mj = markov_algo_parse("""@markovjunior 1 begin
        @calc third() = len / 3
        @grammar begin
            ROOT(third() > 0) => 'R'
            ROOT => 'G'
        end
    end"""),
    result = Ref{String}()
    err = @capture_err (result[] = run_calc(mj))
    @bp_check(result[] == "GGGGGGGGGG", result[])
    @bp_check(n_logged(err) > 0, "Nothing was logged")
end

# Definitions can go anywhere an Op or grammar rule can, and are always global.
let mj = markov_algo_parse("""@markovjunior 1 begin
        @sequence begin
            @calc from_sequence() = 1
            @within_box pixel(min=1, size=1) begin
                @calc from_box() = 2
                @fill 'w'
            end
        end
        @grammar begin
            ROOT => split(x, A(from_sequence()), B(from_box()), C(from_rules()), D(from_leaf()), _)
            @calc from_rules() = 3
            A => 'R'
            B => 'G'
            C => 'B'
            D => begin
                @calc from_leaf() = 4
                @fill 'Y'
            end
        end
    end""")
    @bp_check(collect(keys(mj.calc_macros)) == [ :from_sequence, :from_box, :from_rules, :from_leaf ],
              collect(keys(mj.calc_macros)))
    @bp_check(run_calc(mj) == "RGGBBBYYYY", run_calc(mj))
end
check_parse_fails("""@markovjunior begin
    @calc f() = 1
    @sequence begin
        @calc f() = 2
    end
end""")

# Definitions only check their shape when parsing.
markov_algo_parse("@markovjunior begin\n @calc f() = anything(at_all)\n @fill 'R'\nend")
check_parse_fails("@markovjunior begin\n @calc f() = 1\n @calc f() = 2\nend")
check_parse_fails("@markovjunior begin\n @calc f = 1\nend")
check_parse_fails("@markovjunior begin\n @calc max(a) = a\nend")
check_parse_fails("@markovjunior begin\n @calc\nend")
# Calls in a grammar check their arguments when parsing, but not the function itself.
check_grammar_parse_fails("ROOT => split(x, A(f(nope)))\n A => 'R'")
markov_algo_parse("@markovjunior begin\n @grammar begin\n ROOT => split(x, A(f()))\n A => 'R'\n end\nend")

end)()
