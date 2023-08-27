# Trigonometric functions
import Base.sin, Base.cos, Base.tan, Base.csc, Base.sec, Base.cot

sin(a::GQ) = dimensionless_forward(a, sin)
cos(a::GQ) = dimensionless_forward(a, cos)
tan(a::GQ) = dimensionless_forward(a, tan)

csc(a::GQ) = dimensionless_forward(a, csc)
sec(a::GQ) = dimensionless_forward(a, sec)
cot(a::GQ) = dimensionless_forward(a, cot)

# Inverse trigonometric functions
import Base.asin, Base.acos, Base.atan, Base.acsc, Base.asec, Base.acot

asin(a::GQ) = dimensionless_forward(a, asin)
acos(a::GQ) = dimensionless_forward(a, acos)
atan(a::GQ) = dimensionless_forward(a, atan)

atan(a::GQ, b::GQ) = (a.d == b.d) ? atan(a.x, b.x) : dimension_mismatch((a, b))

acsc(a::GQ) = dimensionless_forward(a, acsc)
asec(a::GQ) = dimensionless_forward(a, asec)
acot(a::GQ) = dimensionless_forward(a, acot)
