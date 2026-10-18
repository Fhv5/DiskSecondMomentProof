module IntervalSafety

using ArbNumerics
using LinearAlgebra

export finite_interval, zero_excluded, verified_inverse, interval_record,
       restore_interval, enclosing_rational, enclosing_rational_interval,
       interval_min, interval_max, rational_bound

"""Require finite endpoints before an interval may justify a decision."""
finite_interval(value::ArbReal) =
    isfinite(value) && isfinite(lowerbound(value)) && isfinite(upperbound(value))

zero_excluded(value::ArbReal) = finite_interval(value) &&
    (upperbound(value) < 0 || lowerbound(value) > 0)

# ArbNumerics 1.6.3 omits the fourth argument in its min/max bindings.
# Keep the workaround local to DSMP and pass FLINT slong precision.

"""Enclose the pointwise minimum with FLINT's required precision argument."""
function interval_min(x::ArbReal{P}, y::ArbReal{P}) where {P}
    result = ArbReal{P}()
    GC.@preserve x y result begin
        ccall(ArbNumerics.@libflint(arb_min), Cvoid,
              (Ref{ArbReal{P}}, Ref{ArbReal{P}}, Ref{ArbReal{P}}, Int),
              result, x, y, P)
    end
    return result
end

"""Enclose the pointwise maximum with FLINT's required precision argument."""
function interval_max(x::ArbReal{P}, y::ArbReal{P}) where {P}
    result = ArbReal{P}()
    GC.@preserve x y result begin
        ccall(ArbNumerics.@libflint(arb_max), Cvoid,
              (Ref{ArbReal{P}}, Ref{ArbReal{P}}, Ref{ArbReal{P}}, Int),
              result, x, y, P)
    end
    return result
end

"""Decode an ArbFloat endpoint to its binary rational value."""
function rational_bound(value::ArbFloat{P}) where {P}
    binary = BigFloat(value; precision=P + 64)
    isfinite(binary) || throw(ArgumentError("nonfinite Arb endpoint"))
    ArbFloat{P}(binary) == value ||
        error("Arb endpoint lost information during binary conversion")
    iszero(binary) && return BigInt(0) // BigInt(1)
    shift = max(0, precision(binary) - exponent(binary) + 2)
    integer = setprecision(BigFloat, precision(binary)) do
        scaled = ldexp(binary, shift)
        isinteger(scaled) || error("Arb endpoint could not be decoded as a binary rational")
        BigInt(scaled)
    end
    return integer // (BigInt(2)^shift)
end

"""Enclose one rational using directed MPFR rounding and verify its bounds."""
function enclosing_rational(::Type{ArbReal{P}}, value::Rational) where {P}
    rational = Rational{BigInt}(value)
    down, up = if rational < 0
        setprecision(BigFloat, P + 2) do
            magnitude = -rational
            (-BigFloat(magnitude, RoundUp; precision=P + 2),
             -BigFloat(magnitude, RoundDown; precision=P + 2))
        end
    else
        (BigFloat(rational, RoundDown; precision=P + 2),
         BigFloat(rational, RoundUp; precision=P + 2))
    end
    down_ball = ArbReal{P}(down)
    up_ball = ArbReal{P}(up)
    result = setinterval(lowerbound(down_ball, ArbFloat),
                         upperbound(up_ball, ArbFloat))
    finite_interval(result) || throw(ArgumentError("nonfinite rational enclosure"))
    rational_bound(lowerbound(result, ArbFloat)) <= rational <=
        rational_bound(upperbound(result, ArbFloat)) ||
        error("Arb interval does not enclose the input rational")
    return result
end

"""Enclose a closed rational interval without rounding either endpoint inward."""
function enclosing_rational_interval(::Type{T}, lower::Rational,
                                     upper::Rational) where {T<:ArbReal}
    lower <= upper || throw(ArgumentError("reversed rational interval"))
    lower_ball = enclosing_rational(T, lower)
    upper_ball = enclosing_rational(T, upper)
    result = setinterval(lowerbound(lower_ball, ArbFloat),
                         upperbound(upper_ball, ArbFloat))
    finite_interval(result) || throw(ArgumentError("nonfinite rational interval"))
    rational_bound(lowerbound(result, ArbFloat)) <= lower &&
        upper <= rational_bound(upperbound(result, ArbFloat)) ||
        error("Arb interval does not enclose the input rational box")
    return result
end

"""Require Arb's determinant enclosure and successful interval-matrix inversion."""
function verified_inverse(matrix::AbstractMatrix{T}) where {T<:ArbReal}
    n, columns = size(matrix)
    n == columns || return nothing
    all(finite_interval, matrix) || return nothing
    work = Matrix{T}(matrix)
    determinant = det(work)
    zero_excluded(determinant) || return nothing
    inverse = try
        inv(work)
    catch err
        err isa ErrorException || rethrow()
        return nothing
    end
    all(finite_interval, inverse) || return nothing
    return inverse, determinant
end

function decimal_rational(value::AbstractString)
    match_result = match(r"^([+-]?)([0-9]+)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$", value)
    match_result === nothing && throw(ArgumentError("unsupported Arb endpoint: $value"))
    fractional = something(match_result.captures[3], "")
    exponent = parse(Int, something(match_result.captures[4], "0")) - length(fractional)
    coefficient = parse(BigInt, match_result.captures[2] * fractional)
    match_result.captures[1] == "-" && (coefficient = -coefficient)
    rational = exponent >= 0 ? coefficient * big(10)^exponent // 1 :
        coefficient // big(10)^(-exponent)
    unit = exponent >= 0 ? big(10)^exponent // 1 : 1 // big(10)^(-exponent)
    return rational, unit
end

"""Encode rational endpoints enlarged until the restored box encloses the original."""
function interval_record(value::ArbReal{P}) where {P}
    finite_interval(value) || throw(ArgumentError("nonfinite interval cannot be archived"))
    approximate_lower, lower_unit = decimal_rational(string(lowerbound(value)))
    approximate_upper, upper_unit = decimal_rational(string(upperbound(value)))
    for shift in 0:2048
        lo = approximate_lower - big(2)^shift * lower_unit
        hi = approximate_upper + big(2)^shift * upper_unit
        record = [string(numerator(lo)), string(denominator(lo)),
                  string(numerator(hi)), string(denominator(hi))]
        restored = restore_interval(ArbReal{P}, record)
        if lowerbound(restored) <= lowerbound(value) &&
           upperbound(value) <= upperbound(restored)
            return record
        end
    end
    error("interval archive failed its outward-containment check")
end

"""Read rational endpoints, then round the reconstructed interval outward."""
function restore_interval(::Type{T}, record) where {T<:ArbReal}
    length(record) == 4 || throw(ArgumentError("expected four rational fields"))
    numbers = parse.(BigInt, String.(record))
    numbers[2] > 0 && numbers[4] > 0 ||
        throw(ArgumentError("interval denominators must be positive"))
    lower = numbers[1] // numbers[2]
    upper = numbers[3] // numbers[4]
    lower <= upper || throw(ArgumentError("reversed interval endpoints"))
    return enclosing_rational_interval(T, lower, upper)
end

# Public API alias.
const exact_rational_bound = rational_bound

end
