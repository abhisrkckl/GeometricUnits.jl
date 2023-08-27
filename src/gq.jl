export GQ

struct GQ{X<:AbstractFloat}
    x::X
    d::Int
end

# Util functions
export dimension_mismatch, dimensionless_forward

dimension_mismatch(obj) = throw(DomainError(obj, "Dimension mismatch!"))
dimensionless_forward(a::GQ, f) = (a.d == 0) ? GQ(f(a.x), 0) : dimension_mismatch(a)

# Unit forwarding
import Base.zero
export unit, quantity_like

unit(a::GQ) = GQ(oneunit(a.x), a.d)
zero(a::GQ) = GQ(zero(a.x), a.d)
quantity_like(a::GQ, y) = GQ(oftype(a.x, y), a.d)

# == Commonly used quantities ===
import Base.time
export dimensionless, speed, distance, frequency, acceleration

dimensionless(x) = GQ(x, 0)
speed(x) = GQ(x, 0)
time(x) = GQ(x, 1)
distance(x) = GQ(x, 1)
frequency(x) = GQ(x, -1)
acceleration(x) = GQ(x, -1)
