# == Addition, subtraction, and negation ===
import Base.+, Base.-

a::GQ + b::GQ = (a.d == b.d) ? GQ(a.x + b.x, a.d) : dimension_mismatch((a, b))
a::GQ + y::Real = (a.d == 0) ? GQ(a.x + y, a.d) : dimension_mismatch((a, y))
x::Real + b::GQ = b + x

a::GQ - b::GQ = (a.d == b.d) ? GQ(a.x - b.x, a.d) : dimension_mismatch((a, b))
a::GQ - y::Real = (a.d == 0) ? GQ(a.x - y, a.d) : dimension_mismatch((a, y))
x::Real - b::GQ = (b.d == 0) ? GQ(x - b.x, b.d) : dimension_mismatch((x, b))

-a::GQ = GQ(-a.x, a.d)

# == Multiplication and division ===
import Base.*, Base./

a::GQ * b::GQ = GQ(a.x * b.x, a.d + b.d)
a::GQ * y::Real = GQ(a.x * y, a.d)
x::Real * b::GQ = b * x

a::GQ / b::GQ = GQ(a.x / b.x, a.d - b.d)
a::GQ / y::Real = GQ(a.x / y, a.d)
x::Real / b::GQ = GQ(x / b.x, -b.d)

# == Power and root ===
import Base.^, Base.inv, Base.sqrt, Base.cbrt
export root

a::GQ^b::GQ = (a.d == 0 && b.d == 0) ? GQ(a.x^b.x, 0) : dimension_mismatch((a, b))
x::Real^b::GQ = (b.d == 0) ? GQ(x^b.x, 0) : dimension_mismatch((x, b))
a::GQ^y::Real = (a.d == 0) ? GQ(a.x^y, 0) : dimension_mismatch((a, y))
a::GQ^n::Integer = GQ(a.x^n, a.d * n)
a::GQ^q::Rational =
    ((a.d * q).den == 1) ? GQ(a.x^q, (a.d * q).num) : dimension_mismatch((a, q))

inv(a::GQ) = GQ(1 / a.x, -a.d)

sqrt(a::GQ) = (a.d % 2 == 0) ? GQ(sqrt(a.x), Int(a.d // 2)) : dimension_mismatch(a)
cbrt(a::GQ) = (a.d % 3 == 0) ? GQ(cbrt(a.x), Int(a.d // 3)) : dimension_mismatch(a)
root(a::GQ, n::Integer) = a^(1 // n)
root(a::GQ, q::Rational) = a^(1 // q)

# == Absolute value ===
import Base.abs

abs(a::GQ) = GQ(abs(a.x), a.d)