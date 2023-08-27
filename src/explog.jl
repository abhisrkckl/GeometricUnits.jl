# Exponential and logarithm
import Base.exp, Base.log, Base.log10, Base.log2

exp(a::GQ) = dimensionless_forward(a, exp)
log(a::GQ) = dimensionless_forward(a, log)
log10(a::GQ) = dimensionless_forward(a, log10)
log2(a::GQ) = dimensionless_forward(a, log2)
