const A = [
    1 2 3
    1 1 1 ;;;

    2 2 2
    1 1 1 ;;;

    1 2 3
    1 1 1
]
const B = [
    1 1 1
    2 2 2
    3 3 3  ;;;

    1 1 1
    2 2 2
    3 3 3  ;;;

    1 1 1
    2 2 2
    3 3 3
]
@bp_check(!MJ.is_symmetric_along(A, 1))
@bp_check(!MJ.is_symmetric_along(A, 2))
@bp_check(MJ.is_symmetric_along(A, 3))
@bp_check(MJ.is_symmetric_along(A, 4))
@bp_check(MJ.is_symmetric_along(A, 5))
@bp_check(!MJ.is_symmetric_along(B, 1))
@bp_check(MJ.is_symmetric_along(B, 2))
@bp_check(MJ.is_symmetric_along(B, 3))

try_swap(a, i, j, flipped, expected) = let result = MJ.are_axes_swappable(a, (i, j), flipped)
    @bp_check(result == expected,
              "Expected are_axes_swappable(", i, ", ", j, ", ", flipped, ") to be ",
                expected, ", but it was ", result)
end
try_swap(A, 1, 2, false, false)
try_swap(A, 1, 3, false, false)
try_swap(A, 1, 5, false, false)
try_swap(A, 6, 6, false, true)
try_swap(A, 1, 2, true, false)
try_swap(A, 1, 3, true, false)
try_swap(A, 1, 5, true, false)
try_swap(A, 6, 6, true, true)
try_swap(B, 1, 2, false, false)
try_swap(B, 1, 2, true, false)
try_swap(B, 2, 3, false, true)
try_swap(B, 2, 3, true, true)
try_swap(B, 1, 3, false, false)
try_swap(B, 1, 3, true, false)
try_swap(B, 1, 6, false, false)
try_swap(B, 1, 6, true, false)
try_swap(B, 1, 6, false, false)
try_swap(B, 7, 6, false, true)
try_swap(B, 7, 6, true, true)