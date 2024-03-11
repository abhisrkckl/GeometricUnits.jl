using ChainRulesCore
using Zygote

dimensionless_zero(a::GQ) = dimensionless(zero(value(a)))
dimensionless_unit(a::GQ) = dimensionless(oneunit(value(a)))

function ChainRulesCore.rrule(::Type{GQ{X}}, x::X, d::Int) where {X<:AbstractFloat}
    y = GQ(x, d)
    function _pullback(ybar)
        fbar = NoTangent()
        xbar = oneunit(y) * ybar
        dbar = NoTangent()
        return fbar, xbar, dbar
    end
    return y, _pullback
end

(p::ProjectTo{X})(a::GQ) where {X<:AbstractFloat} = X(a.x)

Zygote._gradcopy!(dst::AbstractArray, src::GQ) = copyto!(dst, src)

function ChainRulesCore.rrule(::typeof(oneunit), a::GQ)
    y = oneunit(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = dimensionless_zero(a) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(zero), a::GQ)
    y = zero(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = dimensionless_zero(a) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(quantity_like), a::GQ, x)
    y = quantity_like(a, x)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = dimensionless_zero(a) * ybar
        xbar = GQ(oneunit(value(y)), a.d) * ybar
        return fbar, abar, xbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(value), a::GQ)
    y = value(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = 1 / oneunit(a) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(dimensionless), x)
    y = dimensionless(x)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = oneunit(y) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(time), x)
    y = time(x)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = oneunit(y) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(frequency), x)
    y = frequency(x)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = oneunit(x) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(+), a::GQ, b::GQ)
    y = a + b
    function _pullback(ybar)
        fbar = NoTangent()
        abar = dimensionless_unit(a) * ybar
        bbar = dimensionless_unit(b) * ybar
        return fbar, abar, bbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(-), a::GQ, b::GQ)
    y = a - b
    function _pullback(ybar)
        fbar = NoTangent()
        abar = dimensionless_unit(a) * ybar
        bbar = -dimensionless_unit(b) * ybar
        return fbar, abar, bbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(-), a::GQ)
    y = -a
    function _pullback(ybar)
        fbar = NoTangent()
        abar = -dimensionless_unit(a) * ybar
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(*), a::GQ, b::GQ)
    y = a * b
    function _pullback(ybar)
        fbar = NoTangent()
        abar = b * ybar
        bbar = a * ybar
        return fbar, abar, bbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(/), a::GQ, b::GQ)
    y = a / b
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar / b
        bbar = -ybar * a / b / b
        return fbar, abar, bbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(^), a::GQ, b::GQ)
    y = a^b
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * b * y / a
        bbar = ybar * y * log(a)
        return fbar, abar, bbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(^), a::GQ, n::Integer)
    y = a^n
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * signed(n) * y / a
        nbar = NoTangent()
        return fbar, abar, nbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(^), a::GQ, n::Rational)
    y = a^n
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * n * y / a
        nbar = NoTangent()
        return fbar, abar, nbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(sqrt), a::GQ)
    y = sqrt(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * (1 // 2) * y / a
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(cbrt), a::GQ)
    y = cbrt(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * (1 // 3) * y / a
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(root), a::GQ, n::Integer)
    y = root(a, n)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * (1 // n) * y / a
        nbar = NoTangent()
        return fbar, abar, nbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(root), a::GQ, n::Rational)
    y = root(a, n)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * (1 // n) * y / a
        nbar = NoTangent()
        return fbar, abar, nbar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(exp), a::GQ)
    y = exp(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * y
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(log), a::GQ)
    y = log(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar / a
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(log2), a::GQ)
    y = log2(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar / a / log(2)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(log10), a::GQ)
    y = log10(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar / a / log(10)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(sin), a::GQ)
    y = sin(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * cos(a)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(cos), a::GQ)
    y = cos(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = -ybar * sin(a)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(tan), a::GQ)
    y = tan(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * sec(a)^2
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(csc), a::GQ)
    y = csc(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = -ybar * cot(a) * y
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(sec), a::GQ)
    y = sec(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * tan(a) * y
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(cot), a::GQ)
    y = cot(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = -ybar * csc(a)^2
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(asin), a::GQ)
    y = asin(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar / sqrt(1 - a * a)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(acos), a::GQ)
    y = acos(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = -ybar / sqrt(1 - a * a)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(atan), a::GQ)
    y = atan(a)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar / (1 + a * a)
        return fbar, abar
    end
    return y, _pullback
end

function ChainRulesCore.rrule(::typeof(atan), a::GQ, b::GQ)
    y = atan(a, b)
    function _pullback(ybar)
        fbar = NoTangent()
        abar = ybar * b / (a * a + b * b)
        bbar = -ybar * a / (a * a + b * b)
        return fbar, abar, bbar
    end
    return y, _pullback
end
