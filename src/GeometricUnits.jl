module GeometricUnits

import Base.+

dimension_mismatch(obj) = throw(DomainError(obj, "Dimension mismatch!")) 

struct GQ{X <: Real}
    x::X
    d::Int
end 

unit(a::GQ) = typeof(a)(1, a.d)

zero(a::GQ) = typeof(a)(0, a.d)

a::GQ + b::GQ = (a.d == b.d) ? typeof(a)(a.x+b.x, a.d) : dimension_mismatch(obj)

end # module GeometricUnits
