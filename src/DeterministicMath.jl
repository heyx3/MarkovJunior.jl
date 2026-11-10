"
Replacements for important math functions which aren't normally bit-identical on every machine.
The contents of this file allow MarkovJunior.jl to be deterministic across every machine,
  though not necessarily across every Julia version when it comes to random stuff.

(maintained mostly by Claude because I didn't want to rewrite a math standard library myself)

Julia's own `exp`, `log`, `^` etc. use `muladd`, which is fused into a single instruction
  only when the CPU (and the compile target) supports FMA, so their last bits vary between machines.
That's enough to change a whole procedural-generation run, so the engine uses these instead.

They aren't correctly rounded (they're accurate to a few ulps for reasonable inputs);
  what matters is that they're *reproducible*.
They are built only from operations that IEEE 754 defines exactly:
  `+ - * /`, `sqrt`, comparisons, and exact exponent manipulation (`ldexp`, `frexp`, `round`).
Julia never fuses or reorders those on its own (Base's `TwicePrecision` ranges depend on it),
  and `check_determinism()` verifies that while precompiling and on every load.

Rules for editing this file:
* No `muladd`, `fma`, `evalpoly`, `@fastmath`, `@simd`, or Base's transcendental functions.
* Polynomials are written out with Horner's rule, in a fixed order.
* Changing any formula changes results; update `SELF_CHECK_EXPECTED` deliberately when you do.
"
module DeterministicMath

export det_exp, det_exp2, det_log, det_log2, det_log10, det_pow,
       sequential_sum, sequential_prod, sequential_mean, sequential_std


###############
#  Constants

const LN2 = 0.6931471805599453
# `ln(2)` split so that `n * LN2_HI` is exact for any exponent `n` we'll see (from fdlibm).
const LN2_HI = 6.93147180369123816490e-01
const LN2_LO = 1.90821492927058770002e-10
const INV_LN2 = 1.4426950408889634
const INV_LN10 = 0.4342944819032518
const SQRT_HALF = 0.7071067811865476

# Taylor series of `e^r`: 1/k! for k = 0 to 13, enough for |r| <= ln(2)/2.
const EXP_COEFFS = ntuple(k -> 1.0 / Float64(factorial(k - 1)), 14)
# Series of `atanh(s) / s`: 1/(2k+1) for k = 0 to 11, enough for |s| <= 0.1716.
const LOG_COEFFS = ntuple(k -> 1.0 / Float64(2k - 1), 12)


##############
#  Helpers

"`e^r` for small `r` (|r| <= ln(2)/2), by Horner's rule"
function exp_core(r::Float64)::Float64
    p = EXP_COEFFS[14]
    for k in 13:-1:1
        p = (p * r) + EXP_COEFFS[k]
    end
    return p
end

"
Splits a positive finite `x` into `(m, e)` with `x == m * 2^e` and `m` in `[sqrt(1/2), sqrt(2))`,
  then returns `(ln(m), e)`.
"
function log_core(x::Float64)::Tuple{Float64, Float64}
    (m, e) = frexp(x) # m in [0.5, 1)
    if m < SQRT_HALF
        m *= 2.0 # Exact
        e -= 1
    end
    # ln(m) = 2 * atanh(s), where s = (m - 1) / (m + 1).
    s = (m - 1.0) / (m + 1.0)
    z = s * s
    p = LOG_COEFFS[12]
    for k in 11:-1:1
        p = (p * z) + LOG_COEFFS[k]
    end
    return ((2.0 * s) * p, Float64(e))
end

"Raises `x` to a non-negative integer power by repeated squaring, always in the same order"
function pow_by_squaring(x::Float64, n::Int)::Float64
    result = 1.0
    base = x
    while n > 0
        if isodd(n)
            result *= base
        end
        n >>= 1
        if n > 0
            base *= base
        end
    end
    return result
end

is_integer_float(y::Float64) = isfinite(y) && (y == round(y))
is_odd_integer_float(y::Float64) = is_integer_float(y) && (abs(y) < 2.0^53) && isodd(Int64(y))


###############
#  Functions

"`e^x`, reproducible on every machine"
function det_exp(x::Float64)::Float64
    if isnan(x)
        return x
    elseif x > 709.782712893384
        return Inf
    elseif x < -745.1332191019412
        return 0.0
    end
    n = round(x * INV_LN2)
    r = (x - (n * LN2_HI)) - (n * LN2_LO)
    return ldexp(exp_core(r), Int(n))
end

"`2^x`, reproducible on every machine"
function det_exp2(x::Float64)::Float64
    if isnan(x)
        return x
    elseif x >= 1024.0
        return Inf
    elseif x < -1075.0
        return 0.0
    end
    n = round(x)
    f = x - n # Exact
    return ldexp(exp_core(f * LN2), Int(n))
end

"Natural logarithm, reproducible on every machine. Throws a `DomainError` for negative inputs, like Base"
function det_log(x::Float64)::Float64
    if isnan(x) || (x == Inf)
        return x
    elseif x < 0.0
        throw(DomainError(x, "log was called with a negative argument"))
    elseif x == 0.0
        return -Inf
    end
    (ln_m, e) = log_core(x)
    return (e * LN2_HI) + (ln_m + (e * LN2_LO))
end

"Base-2 logarithm, reproducible on every machine. Throws a `DomainError` for negative inputs, like Base"
function det_log2(x::Float64)::Float64
    if isnan(x) || (x == Inf)
        return x
    elseif x < 0.0
        throw(DomainError(x, "log2 was called with a negative argument"))
    elseif x == 0.0
        return -Inf
    end
    (ln_m, e) = log_core(x)
    return e + (ln_m * INV_LN2)
end

"Base-10 logarithm, reproducible on every machine. Throws a `DomainError` for negative inputs, like Base"
det_log10(x::Float64)::Float64 = det_log(x) * INV_LN10

"
`x^y`, reproducible on every machine, following Base's special cases.
Integer powers up to 1024 use repeated squaring (in a fixed order); others use `2^(y * log2(x))`.
Throws a `DomainError` for a negative `x` with a non-integer `y`, like Base.
"
function det_pow(x::Float64, y::Float64)::Float64
    if (y == 0.0) || (x == 1.0)
        return 1.0
    elseif isnan(x) || isnan(y)
        return NaN
    elseif is_integer_float(y) && (abs(y) <= 1024.0)
        n = Int(y)
        return (n >= 0) ? pow_by_squaring(x, n) : (1.0 / pow_by_squaring(x, -n))
    elseif x < 0.0
        if !is_integer_float(y)
            throw(DomainError(x, "Exponentiation yielding a complex result requires a complex argument"))
        end
        magnitude = det_exp2(y * det_log2(-x))
        return is_odd_integer_float(y) ? -magnitude : magnitude
    elseif x == 0.0
        return (y > 0.0) ? 0.0 : Inf
    else
        return det_exp2(y * det_log2(x))
    end
end

# Other number types go through Float64, which is exact for every Float32 and for Ints up to 2^53;
#   rounding the result back to Float32 is exact-rounded too.
for f in (:det_exp, :det_exp2, :det_log, :det_log2, :det_log10)
    @eval $f(x::Real)::Float64 = $f(convert(Float64, x))
    @eval $f(x::Float32)::Float32 = Float32($f(convert(Float64, x)))
end
det_pow(x::Real, y::Real)::Float64 = det_pow(convert(Float64, x), convert(Float64, y))
det_pow(x::Float32, y::Float32)::Float32 = Float32(det_pow(convert(Float64, x), convert(Float64, y)))
"Integer powers are exact (or overflow) in Base's own implementation"
det_pow(x::Integer, y::Integer) = x ^ y


################
#  Reductions

# Base's `sum` and Statistics' `mean`/`std` may add in a different order depending on the CPU's vector width,
#   which changes the last bits of a float result; these always go from left to right.
sequential_sum(a) = foldl(+, a)
sequential_prod(a) = foldl(*, a)
sequential_mean(a) = sequential_sum(a) / length(a)
"The corrected (n-1) standard deviation, like Statistics' `std`"
function sequential_std(a)
    m = sequential_mean(a)
    return sqrt(foldl((total, v) -> begin
        d = v - m
        total + (d * d)
    end, a; init = zero(m)) / (length(a) - 1))
end


################
#  Self-check

# `@noinline` and inputs read from a mutable array, so nothing here is constant-folded:
#   the checks have to run the machine code that's actually used.
@noinline multiply_add(a::Float64, b::Float64, c::Float64) = (a * b) + c
@noinline add_then_subtract(a::Float64, b::Float64) = (a + b) - a

const SELF_CHECK_INPUTS = Float64[
    1.0 + 2.0^-30, 1.0e16,
    1.0, -3.7, 0.3, -10.25,
    10.0, 0.001, 3.0, 7.5,
    0.37, 2.5, 2.0, 0.5, 10.0, -3.3, 0.123, 3.0, -2.0, 5.0
]
"
The exact bits of each `(function, input indices)` case, recorded from this implementation.
If any differ, this machine (or Julia version) would generate different results from everyone else.
"
const SELF_CHECK_CASES = Tuple{Any, Vector{Int}}[
    (det_exp, [ 3 ]), (det_exp, [ 4 ]),
    (det_exp2, [ 5 ]), (det_exp2, [ 6 ]),
    (det_log, [ 7 ]), (det_log, [ 8 ]),
    (det_log2, [ 9 ]), (det_log10, [ 10 ]),
    (det_pow, [ 11, 12 ]), (det_pow, [ 13, 14 ]), (det_pow, [ 15, 16 ]),
    (det_pow, [ 17, 18 ]), (det_pow, [ 19, 20 ])
]
const SELF_CHECK_EXPECTED = UInt64[
    0x4005bf0a8b14576a, 0x3f99511fc6871044,
    0x3ff3b2c47bff8329, 0x3f4ae89f995ad3ae,
    0x40026bb1bbb55516, 0xc01ba18a998fffa0,
    0x3ff95c01a39fbd69, 0x3fec00807a8870ff,
    0x3fb551616d6a40d4, 0x3ff6a09e667f3bcc, 0x3f406c4363887514,
    0x3f5e7d0aba0fa836, 0xc040000000000000
]

"Runs every self-check case, returning the bits of each result"
self_check_results()::Vector{UInt64} = map(SELF_CHECK_CASES) do (f, indices)
    reinterpret(UInt64, f((SELF_CHECK_INPUTS[i] for i in indices)...)::Float64)
end

"
Throws an error if this machine's floating-point math isn't reproducible:
  if Julia fused or reordered plain arithmetic, or any function here gives different bits than recorded.
"
function check_determinism()
    x = SELF_CHECK_INPUTS[1]
    if multiply_add(x, x, -1.0) != 2.0^-29
        error("DeterministicMath: `a * b + c` was fused into a single operation, ",
              "so results won't be reproducible across machines")
    end
    if add_then_subtract(SELF_CHECK_INPUTS[2], 1.0) != 0.0
        error("DeterministicMath: `(a + b) - a` was reordered, ",
              "so results won't be reproducible across machines")
    end
    results = self_check_results()
    for (i, (actual, expected)) in enumerate(zip(results, SELF_CHECK_EXPECTED))
        if actual != expected
            (f, indices) = SELF_CHECK_CASES[i]
            error("DeterministicMath: `", f, "(", join((SELF_CHECK_INPUTS[j] for j in indices), ", "), ")`",
                  " gave bits ", repr(actual), " instead of ", repr(expected),
                  ", so results won't be reproducible across machines")
        end
    end
    return nothing
end

# Check while precompiling, and again on the machine that loads the code.
check_determinism()
__init__() = check_determinism()

end # module
