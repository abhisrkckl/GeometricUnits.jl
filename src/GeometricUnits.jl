module GeometricUnits

import Base.+, Base.-, Base.*, Base./, Base.^, Base.inv, Base.sqrt, Base.cbrt, Base.==, Base.isapprox, Base.time, Base.zero

export GQ,
    dimension_mismatch,
    dimensionless,
    speed,
    time,
    distance,
    frequency,
    acceleration,
    unit,
    zero,
    quantity,
    root

dimension_mismatch(obj) = throw(DomainError(obj, "Dimension mismatch!"))

struct GQ{X<:AbstractFloat}
    x::X
    d::Int
end

# == Unit forwarding ===
unit(a::GQ) = GQ(oneunit(a.x), a.d)
zero(a::GQ) = GQ(zero(a.x), a.d)
quantity(a::GQ, y) = GQ(oftype(a.x, y), a.d)

# == Commonly used quantities ===
dimensionless(x) = GQ(x, 0)
speed(x) = GQ(x, 0)
time(x) = GQ(x, 1)
distance(x) = GQ(x, 1)
frequency(x) = GQ(x, -1)
acceleration(x) = GQ(x, -1)

# == Comparison ===
a::GQ == b::GQ = (a.d == b.d) ? (a.x == b.x) : dimension_mismatch((a, b))
a::GQ == y::Real = (a.d == 0) ? (a.x == y) : dimension_mismatch((a, y))
x::Real == b::GQ = b == x

isapprox(a::GQ, b::GQ) = (a.d == b.d) ? isapprox(a.x, b.x) : dimension_mismatch((a, b))
isapprox(a::GQ, y::Real) = (a.d == 0) ? isapprox(a.x, y) : dimension_mismatch((a, y))
isapprox(x::Real, b::GQ) = isapprox(b, x)

# == Addition, subtraction, and negation ===
a::GQ + b::GQ = (a.d == b.d) ? GQ(a.x + b.x, a.d) : dimension_mismatch((a, b))
a::GQ + y::Real = (a.d == 0) ? GQ(a.x + y, a.d) : dimension_mismatch((a, y))
x::Real + b::GQ = b + x

a::GQ - b::GQ = (a.d == b.d) ? GQ(a.x - b.x, a.d) : dimension_mismatch((a, b))
a::GQ - y::Real = (a.d == 0) ? GQ(a.x - y, a.d) : dimension_mismatch((a, y))
x::Real - b::GQ = (b.d == 0) ? GQ(x - b.x, b.d) : dimension_mismatch((x, b))

-a::GQ = GQ(-a.x, a.d)

# == Multiplication and division ===
a::GQ * b::GQ = GQ(a.x * b.x, a.d + b.d)
a::GQ * y::Real = GQ(a.x * y, a.d)
x::Real * b::GQ = b * x

a::GQ / b::GQ = GQ(a.x / b.x, a.d - b.d)
a::GQ / y::Real = GQ(a.x / y, a.d)
x::Real / b::GQ = GQ(x / b.x, -b.d)

# == Power and root ===
a::GQ ^ b::GQ = (a.d == 0 && b.d == 0) ? GQ(a.x ^ b.x, 0) : dimension_mismatch((a, b))
x::Real ^ b::GQ = (b.d == 0) ? GQ(x ^ b.x, 0) : dimension_mismatch((x, b))
a::GQ ^ y::Real = (a.d == 0) ? GQ(a.x ^ y, 0) : dimension_mismatch((a, y))
a::GQ ^ n::Integer = GQ(a.x ^ n, a.d * n)
a::GQ ^ q::Rational = ((a.d * q).den == 1) ? GQ(a.x ^ q, (a.d * q).num) : dimension_mismatch((a, q))

inv(a::GQ) = GQ(1 / a.x, -a.d)

sqrt(a::GQ) = (a.d % 2 == 0) ? GQ(sqrt(a.x), Int(a.d//2)) : dimension_mismatch(a)
cbrt(a::GQ) = (a.d % 3 == 0) ? GQ(cbrt(a.x), Int(a.d//3)) : dimension_mismatch(a)
root(a::GQ, n::Integer) = a ^ (1 // n)
root(a::GQ, q::Rational) = a ^ (1 // q)

end # module GeometricUnits
