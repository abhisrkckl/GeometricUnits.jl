export taylor_horner, taylor_horner_integral

function validate_taylor_coeffs(x, cs)
    cs_dim = map(udim, cs)
    x_dim = udim(x)
    @assert allequal(x_dim * p + c_dim for (p, c_dim) in enumerate(cs_dim))
end

function taylor_horner(x::GQ, cs)::GQ
    validate_taylor_coeffs(x, cs)

    result_unit = oneunit(first(cs))

    csv = map(value, cs)
    xv = value(x)

    n = length(csv)
    result = zero(last(csv))

    @inbounds for ii = n:-1:1
        result = result * xv / ii + csv[ii]
    end

    return result * result_unit
end

function taylor_horner_integral(x, cs)
    validate_taylor_coeffs(x, cs)

    result_unit = oneunit(first(cs))

    csv = map(value, cs)
    xv = value(x)

    n = length(csv)
    result = zero(last(csv))

    @inbounds for ii = n:-1:1
        result = result * xv / (ii + 1) + csv[ii]
    end

    return (result * result_unit) * x
end

taylor_horner_integral(x, cs, c0) = taylor_horner_integral(x, cs) + c0

taylor_horner(::GQ, cs::Tuple{X}) where {X<:GQ} = cs[1]
taylor_horner(x::GQ, cs::Tuple{X1,X2}) where {X1<:GQ,X2<:GQ} = cs[1] + x * cs[2]
taylor_horner(x::GQ, cs::Tuple{X1,X2,X3}) where {X1<:GQ,X2<:GQ,X3<:GQ} =
    cs[1] + x * (cs[2] + 0.5 * x * cs[3])
taylor_horner_integral(::GQ, cs::Tuple{X}) where {X<:GQ} = x * cs[1]
taylor_horner_integral(x::GQ, cs::Tuple{X1,X2}) where {X1<:GQ,X2<:GQ} =
    x * (cs[1] + 0.5 * x * cs[2])

# function taylor_horner_derivative(x, cs, j)
#     n = length(cs)
#     result = GQ(zero(cs[n].x), cs[n].d - x.d)

#     @inbounds for ii = n:-1:(1+j)
#         result = result * x / (ii - j) + cs[ii]
#     end

#     return result
# end
