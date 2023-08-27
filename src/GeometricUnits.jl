module GeometricUnits

import Base.+, Base.-, Base.*, Base./, Base.==, Base.isapprox, Base.time, Base.zero

export GQ,
    dimension_mismatch,
    dimensionless,
    angle,
    speed,
    time,
    distance,
    frequency,
    acceleration,
    unit,
    zero,
    quantity

dimension_mismatch(obj) = throw(DomainError(obj, "Dimension mismatch!"))

struct GQ{X<:AbstractFloat}
    x::X
    d::Int
end

unit(a::GQ) = GQ(oneunit(a.x), a.d)
zero(a::GQ) = GQ(zero(a.x), a.d)
quantity(a::GQ, y) = GQ(oftype(a.x, y), a.d)

dimensionless(x) = GQ(x, 0)
angle(x) = GQ(x, 0)
speed(x) = GQ(x, 0)
time(x) = GQ(x, 1)
distance(x) = GQ(x, 1)
frequency(x) = GQ(x, -1)
acceleration(x) = GQ(x, -1)

a::GQ == b::GQ = (a.d == b.d) ? (a.x == b.x) : dimension_mismatch((a, b))
a::GQ == y::Real = (a.d == 0) ? (a.x == y) : dimension_mismatch((a, y))
x::Real == b::GQ = b == a

isapprox(a::GQ, b::GQ) = (a.d == b.d) ? isapprox(a.x, b.x) : dimension_mismatch((a, b))
isapprox(a::GQ, y::Real) = (a.d == 0) ? isapprox(a.x, y) : dimension_mismatch((a, y))
isapprox(x::Real, b::GQ) = isapprox(b, x)

a::GQ + b::GQ = (a.d == b.d) ? GQ(a.x + b.x, a.d) : dimension_mismatch((a, b))
a::GQ + y::Real = (a.d == 0) ? GQ(a.x + y, a.d) : dimension_mismatch((a, y))
x::Real + b::GQ = b + x

a::GQ - b::GQ = (a.d == b.d) ? GQ(a.x - b.x, a.d) : dimension_mismatch((a, b))
a::GQ - y::Real = (a.d == 0) ? GQ(a.x - y, a.d) : dimension_mismatch((a, y))
x::Real - b::GQ = (b.d == 0) ? GQ(x - b.x, b.d) : dimension_mismatch((x, b))

-a::GQ = GQ(-a.x, a.d)

a::GQ * b::GQ = GQ(a.x * b.x, a.d + b.d)
a::GQ * y::Real = GQ(a.x * y, a.d)
x::Real * b::GQ = b * x

a::GQ / b::GQ = GQ(a.x / b.x, a.d - b.d)
a::GQ / y::Real = GQ(a.x / y, a.d)
x::Real / b::GQ = GQ(x / b.x, -b.d)

end # module GeometricUnits
