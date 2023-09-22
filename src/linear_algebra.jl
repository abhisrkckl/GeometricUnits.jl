import LinearAlgebra.dot

function dot(xs::Vector{GQ{X}}, ys::Vector{GQ{X}}) where {X<:AbstractFloat}
    res = zero(xs[1] * ys[1])
    for (x, y) in zip(xs, ys)
        res += x * y
    end
    return res
end

function dot(
    xs::Vector{GQ{X}},
    M::Matrix{GQ{X}},
    ys::Vector{GQ{X}},
) where {X<:AbstractFloat}
    p, q = size(M)
    @assert p == length(xs) && q == length(ys) "Shape mismatch in dot."

    res = zero(xs[1] * M[1, 1] * ys[1])

    for i = 1:p
        for j = 1:q
            @inbounds res += xs[i] * M[i, j] * ys[j]
        end
    end

    return res
end
