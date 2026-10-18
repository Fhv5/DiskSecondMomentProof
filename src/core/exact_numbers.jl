"""A number `a + b*sqrt(3)` with rational coefficients."""
struct Qsqrt3 <: Number
    a::Rational{BigInt}
    b::Rational{BigInt}
end

Qsqrt3(a::Rational, b::Rational=0//1) =
    Qsqrt3(Rational{BigInt}(a), Rational{BigInt}(b))
Qsqrt3(a::Integer) = Qsqrt3(a//1, 0//1)
Base.convert(::Type{Qsqrt3}, x::Integer) = Qsqrt3(x)
Base.convert(::Type{Qsqrt3}, x::Rational) = Qsqrt3(x)
Base.promote_rule(::Type{Qsqrt3}, ::Type{T}) where {T<:Integer} = Qsqrt3
Base.promote_rule(::Type{Qsqrt3}, ::Type{Rational{T}}) where {T<:Integer} = Qsqrt3
Base.zero(::Type{Qsqrt3}) = Qsqrt3(0)
Base.one(::Type{Qsqrt3}) = Qsqrt3(1)
Base.iszero(x::Qsqrt3) = iszero(x.a) && iszero(x.b)
Base.:(==)(x::Qsqrt3, y::Qsqrt3) = x.a == y.a && x.b == y.b
Base.:+(x::Qsqrt3, y::Qsqrt3) = Qsqrt3(x.a + y.a, x.b + y.b)
Base.:-(x::Qsqrt3, y::Qsqrt3) = Qsqrt3(x.a - y.a, x.b - y.b)
Base.:-(x::Qsqrt3) = Qsqrt3(-x.a, -x.b)
Base.:*(x::Qsqrt3, y::Qsqrt3) =
    Qsqrt3(x.a*y.a + 3*x.b*y.b, x.a*y.b + x.b*y.a)
Base.inv(x::Qsqrt3) = begin
    denominator = x.a^2 - 3*x.b^2
    iszero(denominator) && throw(DivideError())
    Qsqrt3(x.a/denominator, -x.b/denominator)
end
Base.:/(x::Qsqrt3, y::Qsqrt3) = x * inv(y)
Base.:/(x::Qsqrt3, y::Integer) = Qsqrt3(x.a/y, x.b/y)
Base.:^(x::Qsqrt3, n::Integer) =
    n < 0 ? inv(x^(-n)) : Base.power_by_squaring(x, n)

function _q3_sign(x::Qsqrt3)
    iszero(x) && return 0
    iszero(x.b) && return sign(x.a)
    iszero(x.a) && return sign(x.b)
    sign(x.a) == sign(x.b) && return sign(x.a)
    comparison = x.a^2 - 3*x.b^2
    iszero(comparison) &&
        error("a nonzero rational multiple of sqrt(3) cannot be rational")
    return comparison > 0 ? sign(x.a) : sign(x.b)
end
