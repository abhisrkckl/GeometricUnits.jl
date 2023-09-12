import Base.length

export TaylorSeries

struct TaylorSeries{X<:AbstractFloat}
    cs_rev::Vector{X}
    c0::GQ{X}
    t0::GQ{X}
end

function TaylorSeries(t0::GQ{X}, c0::GQ{X}, cs::Vector{GQ{X}}) where {X<:AbstractFloat}
    _validate_taylor_series(t0, c0, cs)
    cs_rev = reverse([c.x for c in cs])
    return TaylorSeries(cs_rev, c0, t0)
end

function _validate_taylor_series(
    t0::GQ{X},
    c0::GQ{X},
    cs::Vector{GQ{X}},
) where {X<:AbstractFloat}
    dcs = [c.d for c in cs]
    nc = length(cs)
    dc0 = c0.d
    dt0 = t0.d
    dcs1 = dc0 .- (1:nc) .* dt0
    if dcs != dcs1
        dimension_mismatch((cs, c0, t0))
    end
end

length(th::TaylorSeries) = length(th.cs_rev) + 1

function (th::TaylorSeries)(t::GQ)
    result = zero(t.x)
    t_t0 = (t - th.t0).x
    n = length(th)
    for (ii, c) in enumerate(th.cs_rev)
        result = result * t_t0 / (n - ii + 1) + c
        # print(result)
    end
    return quantity_like(th.c0, th.c0.x + t_t0 * result)
end

# a0 + x*(a1 + x*(a2 + a3*x))
