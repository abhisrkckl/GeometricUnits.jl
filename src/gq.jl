export GQ

"""
    GQ{X<:AbstractFloat}

Struct representing a geometrized quantity.
"""
struct GQ{X<:AbstractFloat}
    x::X
    d::Int
end

# Util functions
export dimension_mismatch, dimensionless_forward

"""
    dimension_mismatch(obj)

Utility function to throw a DomainError upon encountering a dimensionally mismatched expression.
"""
dimension_mismatch(obj) = throw(DomainError(obj, "Dimension mismatch!"))

# Unit forwarding
import Base.zero, Base.oneunit
export unit, quantity_like, value, udim

oneunit(a::GQ) = GQ(oneunit(a.x), a.d)
unit(a::GQ) = oneunit(a)
zero(a::GQ) = GQ(zero(a.x), a.d)
quantity_like(a::GQ, y) = GQ(oftype(a.x, y), a.d)
dimensionless_forward(a::GQ, f) = (a.d == 0) ? GQ(f(a.x), 0) : dimension_mismatch((a, f))

value(a::GQ) = a.x
udim(a::GQ) = a.d

# == Commonly used quantities ===
import Base.time
export dimensionless, speed, distance, frequency, acceleration

dimensionless(x) = GQ(x, 0)
speed(x) = GQ(x, 0)
time(x) = GQ(x, 1)
distance(x) = GQ(x, 1)
frequency(x) = GQ(x, -1)
acceleration(x) = GQ(x, -1)
