import Base.iterate
import LinearAlgebra.dot, LinearAlgebra.transpose

iterate(a::GQ) = (a, nothing)
iterate(::GQ, ::Any) = nothing

transpose(a::GQ) = a

dot(a::GQ, b::GQ) = a * b
