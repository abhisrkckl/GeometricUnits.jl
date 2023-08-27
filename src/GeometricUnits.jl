module GeometricUnits

import Base.+, Base.-, Base.time, Base.zero

export GQ, dimension_mismatch, dimensionless, angle, speed, time, distance, frequency, acceleration, unit, zero, quantity

dimension_mismatch(obj) = throw(DomainError(obj, "Dimension mismatch!"))

struct GQ{X<:Real}
    x::X
    d::Int
end

unit(a::GQ) = GQ(oneunit(a.x), a.d)
zero(a::GQ) = GQ(zero(a.x), a.d)
quantity(a::GQ, y) = GQ(oftype(a.x, y), a.d)

dimensionless(x::Float64) = GQ{Float64}(x, 0)
angle(x::Float64) = GQ{Float64}(x, 0)
speed(x::Float64) = GQ{Float64}(x, 0)

time(x::Float64) = GQ{Float64}(x, 1)
distance(x::Float64) = GQ{Float64}(x, 1)

frequency(x::Float64) = GQ{Float64}(x, -1)
acceleration(x::Float64) = GQ{Float64}(x, -1)

a::GQ + b::GQ = (a.d == b.d) ? GQ(a.x + b.x, a.d) : dimension_mismatch((a, b))

a::GQ - b::GQ = (a.d == b.d) ? GQ(a.x - b.x, a.d) : dimension_mismatch((a, b))

-a::GQ = GQ(-a.x, a.d)

end # module GeometricUnits
