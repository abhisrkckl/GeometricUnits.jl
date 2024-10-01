# Exponential and logarithm
import Base.exp, Base.exp10, Base.exp2, Base.log, Base.log10, Base.log2

exp(a::GQ{0,X}) where {X} = dimensionless(exp(a.x))
exp10(a::GQ{0,X}) where {X} = dimensionless(exp10(a.x))
exp2(a::GQ{0,X}) where {X} = dimensionless(exp2(a.x))
log(a::GQ{0,X}) where {X} = dimensionless(log(a.x))
log10(a::GQ{0,X}) where {X} = dimensionless(log10(a.x))
log2(a::GQ{0,X}) where {X} = dimensionless(log2(a.x))
