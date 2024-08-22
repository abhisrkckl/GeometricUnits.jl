# == Addition, subtraction, and negation ===
import Base.+, Base.-

a::GQ + b::GQ = (a.d == b.d) ? GQ(a.x + b.x, a.d) : dimension_mismatch((a, b))
a::GQ + y::Real = a + dimensionless(float(y))
x::Real + b::GQ = b + x

a::GQ - b::GQ = (a.d == b.d) ? GQ(a.x - b.x, a.d) : dimension_mismatch((a, b))
a::GQ - y::Real = a - dimensionless(float(y))
x::Real - b::GQ = -(b - x)

-a::GQ = GQ(-a.x, a.d)

# == Multiplication and division ===
import Base.*, Base./

a::GQ * b::GQ = GQ(a.x * b.x, a.d + b.d)
a::GQ * y::Real = a * dimensionless(float(y))
x::Real * b::GQ = b * x

a::GQ / b::GQ = GQ(a.x / b.x, a.d - b.d)
a::GQ / y::Real = a / dimensionless(float(y))
x::Real / b::GQ = dimensionless(float(x)) / b

# == Power and root ===
import Base.^, Base.inv, Base.sqrt, Base.cbrt
export root

a::GQ^b::GQ = (a.d == 0 && b.d == 0) ? GQ(a.x^b.x, 0) : dimension_mismatch((a, b))
x::Real^b::GQ = dimensionless(float(x))^b
a::GQ^y::Real = a^dimensionless(float(y))
a::GQ^n::Integer = GQ(a.x^n, a.d * signed(n))
a::GQ^q::Rational =
    ((a.d * q).den == 1) ? GQ(a.x^q, (a.d * q).num) : dimension_mismatch((a, q))

inv(a::GQ) = 1 / a

sqrt(a::GQ) = (a.d % 2 == 0) ? GQ(sqrt(a.x), Int(a.d // 2)) : dimension_mismatch(a)
cbrt(a::GQ) = (a.d % 3 == 0) ? GQ(cbrt(a.x), Int(a.d // 3)) : dimension_mismatch(a)
root(a::GQ, n::Union{Integer,Rational}) = a^(1 // n)

# == Absolute value ===
import Base.abs

abs(a::GQ) = GQ(abs(a.x), a.d)

# == Sign ==============
import Base.sign

sign(a::GQ) = sign(a.x)