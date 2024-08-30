import Base.iterate, Base.length
import LinearAlgebra.dot, LinearAlgebra.transpose

transpose(a::GQ) = a

dot(a::GQ, b::GQ) = a * b
dot(a::GQ, b::Number) = a * b
dot(a::Number, b::GQ) = a * b
