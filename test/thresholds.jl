# Tests for thresholds, which are calculations (see `CalcInterpreter`).

(() -> begin

"Runs the algorithm on a 1D grid of black cells, returning how many became white"
function count_white(mj::MarkovAlgorithm; len = 20, seed = 0x1234)
    # Thresholds always print as `threshold(...)`, and must survive the round-trip.
    @bp_check(markov_algo_parse(markov_algo_to_string(mj)) == mj,
              "Round-trip failed:\n", markov_algo_to_string(mj))
    state_ch = markov_algo_run(mj, fill(MJ.CELL_CODE_BY_CHAR['b'], len),
                               MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                               seeds=seed)
    return markov_algo_complete(g -> count(==(MJ.CELL_CODE_BY_CHAR['w']), g), state_ch)
end
count_white(src::String; kw...) = count_white(markov_algo_parse("@markovjunior 1 begin\n$src\nend"); kw...)

# The classic forms still work.
@bp_check(count_white("@rewrite 3 b => w") == 3)
@bp_check(count_white("@rewrite (size_prod()/4) b => w") == 5)
@bp_check(count_white("@rewrite (area/4) b => w") == 5) # `area` is short for `size_prod()`
@bp_check(count_white("@rewrite (0.5*size_average()) b => w") == 10)
@bp_check(count_white("@rewrite (size_prod()*0.1) b => w") == 2)
@bp_check(count_white("@rewrite b => w") == 20)
@bp_check(count_white("@rewrite (size_prod()/100):20 b => w") in 0:20)
# Mixed Int/Float ranges are Float ranges, and the result is rounded.
let seen = Set(count_white("@rewrite (2:4.0) b => w"; seed = s) for s in 1:60)
    @bp_check(seen == Set(2:4), seen)
end

# Any calculation works, including user functions.
@bp_check(count_white("@rewrite (clamp(size_prod() ÷ 3, 1, 4)) b => w") == 4)
@bp_check(count_white("@rewrite (size(x) / 4) b => w") == 5) # Axes can be named
@bp_check(count_white("@calc wanted() = 7\n@rewrite (wanted()) b => w") == 7)
@bp_check(count_white("@calc wanted() = 7\n@rewrite wanted() b => w") == 7) # No parentheses needed without spaces
@bp_check(count_white("@calc wanted() = 7\n@rewrite threshold(wanted()) b => w") == 7)
# The rules-then-biases form is still recognized.
@bp_check(count_white("@rewrite b => w temperature(0.5)") == 20)
@bp_check(count_white("@rewrite begin\n b => w\n end temperature(0.5)") == 20)

# Results are rounded, and 0 or less does nothing.
@bp_check(count_white("@rewrite (2.6) b => w") == 3)
@bp_check(count_white("@rewrite (0) b => w") == 0)
@bp_check(count_white("@rewrite (-5) b => w") == 0)
@bp_check(count_white("@sequence (0) begin
 @rewrite 1 b => w
 end") == 0)
@bp_check(count_white("@sequence (-1) begin
 @rewrite 1 b => w
 end") == 0)

# Runtime problems are logged, and fall back to 1.
for threshold in ("nope()", "(1 / 0)", "(size_prod() > 3)")
    result = Ref{Int}()
    err = @capture_err (result[] = count_white("@rewrite $threshold b => w"))
    @bp_check(result[] == 1, threshold, " gave ", result[])
    @bp_check(!isempty(err), "Nothing was logged for ", threshold)
end

# User functions can be replaced between runs.
let mj = markov_algo_parse("""@markovjunior 1 begin
        @calc wanted() = 2
        @rewrite wanted() b => w
    end""")
    @bp_check(count_white(mj) == 2)
    push!(mj.calc_macros, calc_macro_parse("wanted() = size_prod() - 1"))
    @bp_check(count_white(mj) == 19)
end

# `@sequence`: the first of two arguments is the threshold unless it's an Op or a block of them.
@bp_check(count_white("@sequence 3 @rewrite 1 b => w") == 3)
@bp_check(count_white("@calc reps() = 4\n@sequence reps() @rewrite 1 b => w") == 4)
@bp_check(count_white("@sequence threshold(2) @rewrite 1 b => w") == 2)
@bp_check(count_white("@sequence (size_average()/10) begin\n @rewrite 1 b => w\n end") == 2)
@bp_check(count_white("@sequence begin\n @rewrite 1 b => w\n end temperature(0.5)") == 1)
@bp_check(count_white("@sequence repeat @rewrite 2 b => w") == 20)

# `@within_regions`: anything that isn't `inherit_biases`, `border(...)`, or `along(...)` is the threshold.
let src = """@fill 'w' pixel(min=5, size=1)
             @fill 'w' pixel(min=10, size=1)
             @fill 'w' pixel(min=15, size=1)
             @calc picks() = 2
             @within_regions border(0) picks() begin
                 b(g) => @fill 'w'
             end"""
    @bp_check(count_white(src) in (3 + 8, 3 + 9), count_white(src)) # 2 of the 4 black gaps (4 cells, 4, 4, 5) get filled
end

end)()
