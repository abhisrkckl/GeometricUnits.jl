# == Comparison ===
import Base.==,
    Base.isapprox, Base.<, Base.<=, Base.>, Base.>=, Base.isfinite, Base.isinf, Base.isnan

(a::GQ{d,X} == b::GQ{d,Y}) where {d,X,Y} = (a.x == b.x)
(a::GQ{0,X} == y::Real) where {X} = (a.x == y)
(x::Real == b::GQ) = (b == x)

isapprox(a::GQ{d,X}, b::GQ{d,Y}) where {d,X,Y} = isapprox(a.x, b.x)
isapprox(a::GQ{0,X}, y::Real) where {X} = isapprox(a.x, y)
isapprox(x::Real, b::GQ) = isapprox(b, x)

(a::GQ{d,X} > b::GQ{d,Y}) where {d,X,Y} = (a.x > b.x)
(a::GQ{0,X} > y::Real) where {X} = (a.x > y)
(x::Real > b::GQ{0,X}) where {X} = (x > b.x)

(a::GQ{d,X} >= b::GQ{d,Y}) where {d,X,Y} = (a.x >= b.x)
(a::GQ{0,X} >= y::Real) where {X} = (a.x >= y)
(x::Real >= b::GQ{0,X}) where {X} = (x >= b.x)

(a::GQ{d,X} < b::GQ{d,Y}) where {d,X,Y} = (a.x < b.x)
(a::GQ{0,X} < y::Real) where {X} = (a.x < y)
(x::Real < b::GQ{0,X}) where {X} = (x < b.x)

(a::GQ{d,X} <= b::GQ{d,Y}) where {d,X,Y} = (a.x <= b.x)
(a::GQ{0,X} <= y::Real) where {X} = (a.x <= y)
(x::Real <= b::GQ{0,X}) where {X} = (x <= b.x)

isfinite(a::GQ) = isfinite(a.x)
isnan(a::GQ) = isnan(a.x)
isinf(a::GQ) = isinf(a.x)
