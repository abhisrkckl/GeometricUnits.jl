# == Addition, subtraction, and negation ===
import Base.+, Base.-

(a::GQ{d,X} + b::GQ{d,Y}) where {d,X,Y} = GQ{d}(a.x + b.x)
(a::GQ{0,X} + y::Real) where {X} = a + dimensionless(float(y))
(x::Real + b::GQ) = (b + x)

(a::GQ{d,X} - b::GQ{d,Y}) where {d,X,Y} = GQ{d}(a.x - b.x)
(a::GQ{0,X} - y::Real) where {X} = a - dimensionless(float(y))
x::Real - b::GQ = -(b - x)

(-a::GQ{d,X}) where {d,X} = GQ{d}(-a.x)
+a::GQ = a

# == Multiplication and division ===
import Base.*, Base./

(a::GQ{d1,X} * b::GQ{d2,Y}) where {d1,d2,X,Y} = GQ{d1+d2}(a.x * b.x)
a::GQ * y::Real = a * dimensionless(float(y))
x::Real * b::GQ = b * x

(a::GQ{d1,X} / b::GQ{d2,Y}) where {d1,d2,X,Y} = GQ{d1-d2}(a.x / b.x)
a::GQ / y::Real = a / dimensionless(float(y))
x::Real / b::GQ = dimensionless(float(x)) / b

# == Power and root ===
import Base.^, Base.inv, Base.sqrt, Base.cbrt
export root

(a::GQ{0,X}^b::GQ{0,Y}) where {X,Y} = dimensionless(a.x^b.x)
(x::Real^b::GQ{0,Y}) where {Y} = dimensionless(float(x))^b
(a::GQ{0,X}^y::Real) where {X} = a^dimensionless(float(y))
(a::GQ{d,X}^vq::Val{q}) where {d,q,X} = GQ{Int(q*d)}(a.x^q)

inv(a::GQ) = 1 / a

sqrt(a::GQ{d,X}) where {d,X} = GQ{Int(d//2)}(sqrt(a.x))
cbrt(a::GQ{d,X}) where {d,X} = GQ{Int(d//3)}(cbrt(a.x))
root(a::GQ{0,X}, q::Real) where {X} = a^(1 // q)
root(a::GQ, ::Val{q}) where {q} = a^Val(1 // q)

# == Absolute value ===
import Base.abs

abs(a::GQ) = GQ(abs(a.x), a.d)

# == Sign ==============
import Base.sign

sign(a::GQ) = sign(a.x)

# == Floor and ceiling =
import Base.floor, Base.ceil

floor(a::GQ) = GQ(floor(a.x), a.d)
ceil(a::GQ) = GQ(ceil(a.x), a.d)