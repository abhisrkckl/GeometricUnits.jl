export GQ

"""
    GQ{d,X<:AbstractFloat}

Represents a quantity with dimensions ``[T^d]``.
A `GQ` can be represented mathematically as ``x s^d`` where ``d ∈ ℤ``.
"""
struct GQ{d,X<:AbstractFloat}
    x::X

    function GQ{d,X}(x::X) where {d,X<:AbstractFloat}
        @assert d isa Int
        return new(x)
    end
end

GQ{d}(x::X) where {d,X<:AbstractFloat} = GQ{d,X}(x)
GQ{X}(a::GQ{d,Y}) where {d,X,Y} = GQ{d}(X(a.x))

# Extract value and dimensions
export value, udim

value(a::GQ) = a.x
udim(::GQ{d,X}) where {d,X} = d

# Unit forwarding and conversion
import Base.zero, Base.oneunit, Base.oftype, Base.convert

oneunit(a::GQ{d,X}) where {d,X} = GQ{d}(oneunit(a.x))
zero(a::GQ{d,X}) where {d,X} = GQ{d}(zero(a.x))
oftype(a::GQ{d,X}, y) where {d,X} = GQ{d}(oftype(a.x, y))
convert(::Type{X}, a::GQ{d,X}) where {d,X} = a.x

oneunit(::Type{GQ{d,X}}) where {d,X} = GQ{d}(oneunit(X))
zero(::Type{GQ{d,X}}) where {d,X} = GQ{d}(zero(X))

# == Iteration ==================
import Base.iterate, Base.length

iterate(a::GQ) = (a, nothing)
iterate(::GQ, ::Any) = nothing
length(::GQ) = 1

# == Commonly used quantities ===
import Base.time
export dimensionless, speed, distance, mass, frequency, acceleration

dimensionless(x) = GQ{0}(x)
speed(x) = dimensionless(x)
time(x) = GQ{1}(x)
distance(x) = time(x)
mass(x) = time(x)
frequency(x) = GQ{-1}(x)
acceleration(x) = frequency(x)
