# == Comparison ===
import Base.==, Base.isapprox, Base.<, Base.<=, Base.>, Base.>=, Base.isfinite, Base.isinf, Base.isnan

a::GQ == b::GQ = (a.d == b.d) ? (a.x == b.x) : dimension_mismatch((a, b))
a::GQ == y::Real = (a.d == 0) ? (a.x == y) : dimension_mismatch((a, y))
x::Real == b::GQ = b == x

isapprox(a::GQ, b::GQ) = (a.d == b.d) ? isapprox(a.x, b.x) : dimension_mismatch((a, b))
isapprox(a::GQ, y::Real) = (a.d == 0) ? isapprox(a.x, y) : dimension_mismatch((a, y))
isapprox(x::Real, b::GQ) = isapprox(b, x)

a::GQ > b::GQ = (a.d == b.d) ? (a.x > b.x) : dimension_mismatch((a, b))
a::GQ > y::Real = (a.d == 0) ? (a.x > y) : dimension_mismatch((a, y))
x::Real > b::GQ = b < x

a::GQ >= b::GQ = (a.d == b.d) ? (a.x >= b.x) : dimension_mismatch((a, b))
a::GQ >= y::Real = (a.d == 0) ? (a.x >= y) : dimension_mismatch((a, y))
x::Real >= b::GQ = b <= x

a::GQ < b::GQ = (a.d == b.d) ? (a.x < b.x) : dimension_mismatch((a, b))
a::GQ < y::Real = (a.d == 0) ? (a.x < y) : dimension_mismatch((a, y))
x::Real < b::GQ = b > x

a::GQ <= b::GQ = (a.d == b.d) ? (a.x <= b.x) : dimension_mismatch((a, b))
a::GQ <= y::Real = (a.d == 0) ? (a.x <= y) : dimension_mismatch((a, y))
x::Real <= b::GQ = b >= x

isfinite(a::GQ) = isfinite(a.x)
isnan(a::GQ) = isnan(a.x)
isinf(a::GQ) = isinf(a.x)