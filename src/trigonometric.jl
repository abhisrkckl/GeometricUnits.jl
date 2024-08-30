# Trigonometric functions
import Base.sin, Base.cos, Base.sincos, Base.tan, Base.csc, Base.sec, Base.cot

sin(a::GQ{0,X}) where {X} = dimensionless(sin(a.x))
cos(a::GQ{0,X}) where {X} = dimensionless(cos(a.x))
sincos(a::GQ{0,X}) where {X} = map(dimensionless, sincos(a.x))
tan(a::GQ{0,X}) where {X} = dimensionless(tan(a.x))

csc(a::GQ{0,X}) where {X} = dimensionless(csc(a.x))
sec(a::GQ{0,X}) where {X} = dimensionless(sec(a.x))
cot(a::GQ{0,X}) where {X} = dimensionless(cot(a.x))

# Inverse trigonometric functions
import Base.asin, Base.acos, Base.atan, Base.acsc, Base.asec, Base.acot

asin(a::GQ{0,X}) where {X} = dimensionless(asin(a.x))
acos(a::GQ{0,X}) where {X} = dimensionless(acos(a.x))
atan(a::GQ{0,X}) where {X} = dimensionless(atan(a.x))

atan(a::GQ{d,X}, b::GQ{d,Y}) where {d,X,Y} = dimensionless(atan(a.x, b.x))

acsc(a::GQ{0,X}) where {X} = dimensionless(acsc(a.x))
asec(a::GQ{0,X}) where {X} = dimensionless(asec(a.x))
acot(a::GQ{0,X}) where {X} = dimensionless(acot(a.x))
