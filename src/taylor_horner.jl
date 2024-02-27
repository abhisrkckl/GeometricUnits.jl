export taylor_horner, taylor_horner_integral

function taylor_horner(x, cs)
    n = length(cs)
    result = GQ(zero(cs[n].x), cs[n].d - x.d)

    @inbounds for ii = n:-1:1
        result = result * x / ii + cs[ii]
    end

    return result
end

function taylor_horner_integral(x, cs, c0)
    n = length(cs)
    result = GQ(zero(c0.x), cs[n].d - x.d)

    @inbounds for ii = n:-1:1
        result = result * x / (ii + 1) + cs[ii]
    end

    return result * x + c0
end
