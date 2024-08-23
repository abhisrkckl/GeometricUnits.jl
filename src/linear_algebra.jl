import Base.iterate, Base.length
import LinearAlgebra.dot, LinearAlgebra.transpose

iterate(a::GQ) = (a, nothing)
iterate(::GQ, ::Any) = nothing

transpose(a::GQ) = a

dot(a::GQ, b::GQ) = a * b
dot(a::GQ, b::Number) = a * b
dot(a::Number, b::GQ) = a * b

length(::GQ) = 1