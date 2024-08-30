# Exponential and logarithm
import Base.exp, Base.log, Base.log10, Base.log2

exp(a::GQ{0,X}) where {X} = dimensionless(exp(a.x))
log(a::GQ{0,X}) where {X} = dimensionless(log(a.x))
log10(a::GQ{0,X}) where {X} = dimensionless(log10(a.x))
log2(a::GQ{0,X}) where {X} = dimensionless(log2(a.x))
