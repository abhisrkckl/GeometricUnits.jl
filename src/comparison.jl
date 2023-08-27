# == Comparison ===
import Base.==, Base.isapprox

a::GQ == b::GQ = (a.d == b.d) ? (a.x == b.x) : dimension_mismatch((a, b))
a::GQ == y::Real = (a.d == 0) ? (a.x == y) : dimension_mismatch((a, y))
x::Real == b::GQ = b == x

isapprox(a::GQ, b::GQ) = (a.d == b.d) ? isapprox(a.x, b.x) : dimension_mismatch((a, b))
isapprox(a::GQ, y::Real) = (a.d == 0) ? isapprox(a.x, y) : dimension_mismatch((a, y))
isapprox(x::Real, b::GQ) = isapprox(b, x)
