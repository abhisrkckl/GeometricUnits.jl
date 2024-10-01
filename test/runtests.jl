using GeometricUnits
using LinearAlgebra
using Quadmath
using Test
using TestExtras
using BenchmarkTools
using Zygote

@testset verbose = true begin
    t1 = time(2.9)
    t2 = oftype(t1, 1.2)

    f1 = frequency(0.2)

    d1 = dimensionless(6.0)
    a1 = dimensionless(0.2)

    @testset "constructor" begin
        t11 = GQ{Float32}(t1)
        @test isa(t11.x, Float32)
        @test value(t11) == Float32(value(t1))
        @test udim(t11) == udim(t1)
        @test length(t11) == 1
        @test_throws AssertionError GQ{1.1}(1.2)
    end

    @testset "unit and zero" begin
        @test value(oneunit(t1)) == 1 && udim(oneunit(t1)) == udim(t1)
        @test value(zero(t1)) == 0 && udim(zero(t1)) == udim(t1)

        @test value(oneunit(f1)) == 1 && udim(oneunit(f1)) == udim(f1)
        @test value(zero(f1)) == 0 && udim(zero(f1)) == udim(f1)

        @constinferred oneunit(t1)
        @constinferred value(t1)
        @constinferred udim(t1)
        @constinferred zero(t1)
    end

    @testset "commonly used quantities" begin
        @test dimensionless(2.5) == speed(2.5)
        @test time(2.5) == distance(2.5)
        @test frequency(2.5) == acceleration(2.5)

        @constinferred dimensionless(2.5)
        @constinferred speed(2.5)
        @constinferred time(2.5)
        @constinferred distance(2.5)
        @constinferred frequency(2.5)
        @constinferred acceleration(2.5)
        @constinferred mass(2.1)
    end

    @testset "comparison" begin
        @test t1 == t1
        @test !(t1 != t1)
        @test t1 != f1

        @constinferred (t1 == t1)
        @constinferred (t1 != t1)

        @test t1 ≈ t1
        @test_throws MethodError t1 ≈ f1

        @constinferred (t1 ≈ t1)

        @test 0 == zero(d1)
        @test 0 ≈ zero(d1)

        @constinferred (0 == zero(d1))
        @constinferred (0 ≈ zero(d1))

        @test (t1 > t2) || (t1 <= t2)
        @test (t1 < t2) || (t1 >= t2)

        @constinferred (t1 > t2)
        @constinferred (t1 <= t2)
        @constinferred (t1 < t2)
        @constinferred (t1 >= t2)

        @test (d1 > 0) || (d1 <= 0)
        @test (d1 < 0) || (d1 >= 0)
        @test (0 > d1) || (0 <= d1)
        @test (0 < d1) || (0 >= d1)

        @constinferred (d1 > 0)
        @constinferred (d1 <= 0)
        @constinferred (d1 < 0)
        @constinferred (d1 >= 0)
        @constinferred (0 > d1)
        @constinferred (0 <= d1)
        @constinferred (0 < d1)
        @constinferred (0 >= d1)

        @test_throws MethodError d1 > t1
        @test_throws MethodError d1 >= t1
        @test_throws MethodError d1 <= t1
        @test_throws MethodError d1 < t1

        @test_throws MethodError 0 > t1
        @test_throws MethodError 0 >= t1
        @test_throws MethodError 0 <= t1
        @test_throws MethodError 0 < t1

        @test isfinite(t1)
        @test !isnan(t1)
        @test !isinf(t1)

        @constinferred isfinite(t1)
        @constinferred isnan(t1)
        @constinferred isinf(t1)
    end

    @testset "addition and subtraction" begin
        @test (t1 + t2).x == t1.x + t2.x
        @test_throws MethodError t1 + f1

        @test (t1 - t2).x == t1.x - t2.x
        @test_throws MethodError t1 - f1

        @test (t1 - t2).x == (t1 + (-t2)).x
        @test (t1 + t2).x == (t1 - (-t2)).x

        @test t1 + t2 == t2 + t1
        @test t1 - t2 == -(t2 - t1)
        @test t1 - t1 == zero(t1)

        @test 1 + d1 == +(d1 + 1)
        @test 1 - d1 == -(d1 - 1)

        @test_throws MethodError 1 + t1
        @test_throws MethodError 1 - t1
        @test_throws MethodError t1 + 1
        @test_throws MethodError t1 - 1

        @constinferred (t1 + t2)
        @constinferred (t1 - t2)
        @constinferred (-t1)
        @constinferred (1 + d1)
        @constinferred (d1 + 1)
        @constinferred (1 - d1)
        @constinferred (d1 - 1)
    end

    @testset "multiplication and division" begin
        @test (t1 * f1).x == t1.x * f1.x
        @test udim(t1 * f1) == udim(t1) + udim(f1)

        @test t1 * d1 == d1 * t1
        @test t1 * d1 == d1.x * t1
        @test t1 * d1 == t1 * d1.x

        @test t1 + t1 == 2 * t1
        @test t1 + 3.2 * t1 == 4.2 * t1

        @test (t1 / f1).x == t1.x / f1.x
        @test udim(t1 / f1) == udim(t1) - udim(f1)

        @test t1 / f1 ≈ t1 * (1 / f1)
        @test t1 / t1 == 1

        @test t1 + 1 / f1 ≈ (t1 * f1 + 1) / f1

        @test 1 / (1 / d1) ≈ d1
        @test 2 * (d1 / 2) ≈ d1

        @constinferred (t1 * f1)
        @constinferred (t1 / f1)
        @constinferred (1.2 * f1)
        @constinferred (f1 * 1.3)
        @constinferred (1.2 / f1)
        @constinferred (f1 / 1.3)
        @constinferred map(*, Ref(d1), (t1, f1))
    end

    @testset "power and root" begin
        @test t1 * t1 == t1^Val(2)
        @test 1 / t1 == t1^Val(-1)

        @test sqrt(t1^Val(2)) ≈ t1
        @test cbrt(t1^Val(3)) ≈ t1
        @test root(t1^Val(5), Val(5)) ≈ t1

        @test sqrt(d1) ≈ d1^Val(0.5)
        @test sqrt(d1) ≈ d1^Val(1 // 2)
        @test d1^Val(2) ≈ root(d1, Val(1 // 2))

        @test 1.0^d1 ≈ 1
        @test d1^0 ≈ 1

        @test d1^d1 == d1.x^d1.x

        @test_throws MethodError t1^(1 // 2)
        @test_throws MethodError t1^6.3
        @test_throws MethodError 0.5^t1
        @test_throws InexactError sqrt(t1)
        @test_throws InexactError cbrt(t1)
        @test_throws MethodError root(t1, 4)

        @constinferred (t1^Val(2))
        @constinferred (2^d1)
        @constinferred sqrt(d1)
        @constinferred cbrt(d1)
        @constinferred root(d1, 4)
    end

    @testset "abs, sign, ceil, and floor" begin
        @test abs(d1) == d1
        @test abs(-d1) == d1

        @test sign(d1) == 1
        @test sign(-d1) == -1
        @test sign(zero(d1)) == 0

        @test floor(frequency(1.1)) == frequency(1.0)
        @test ceil(frequency(1.1)) == frequency(2.0)
    end

    @testset "exp and log" begin
        @test exp(log(d1)) ≈ d1
        @test exp10(log10(d1)) ≈ d1
        @test exp2(log2(d1)) ≈ d1

        @constinferred exp(d1)
        @constinferred exp10(d1)
        @constinferred exp2(d1)
        @constinferred log(d1)
        @constinferred log10(d1)
        @constinferred log2(d1)
    end

    @testset "trigonometric functions" begin
        @test sin(d1) ≈ 1 / csc(d1)
        @test cos(d1) ≈ 1 / sec(d1)
        @test tan(d1) ≈ 1 / cot(d1)
        @test sincos(d1) == (sin(d1), cos(d1))

        @test asin(sin(a1)) ≈ a1
        @test acos(cos(a1)) ≈ a1
        @test atan(tan(a1)) ≈ a1
        @test acsc(csc(a1)) ≈ a1
        @test asec(sec(a1)) ≈ a1
        @test acot(cot(a1)) ≈ a1

        @test atan(t1, t2) ≈ atan(t1 / t2)

        @constinferred sin(d1)
        @constinferred csc(d1)
        @constinferred cos(d1)
        @constinferred sec(d1)
        @constinferred tan(d1)
        @constinferred cot(d1)

        @constinferred asin(sin(a1))
        @constinferred acos(cos(a1))
        @constinferred atan(tan(a1))
        @constinferred acsc(csc(a1))
        @constinferred asec(sec(a1))
        @constinferred acot(cot(a1))

        @constinferred atan(t1, t2)
        @constinferred atan(d1)
    end

    @testset "taylor-horner" begin
        t0 = time(0.0)
        c0 = distance(1.5)
        c1 = speed(1.0)
        c2 = acceleration(1.2)
        c3 = acceleration(0.5) / time(2.0)
        c4 = acceleration(1.2) / time(2.0)^Val(2)
        cs = (c1, c2, c3, c4)

        t1 = time(1.0)

        dt = t1 - t0

        @test taylor_horner(dt, cs) ≈
              (c1 + c2 * dt + (1 / 2) * c3 * dt * dt + (1 / 6) * c4 * dt * dt * dt)

        @test (@ballocated taylor_horner($dt, $cs)) == 0

        @test taylor_horner_integral(dt, cs, c0) ≈ (
            c0 +
            c1 * dt +
            (1 / 2) * c2 * dt^Val(2) +
            (1 / 6) * c3 * dt^Val(3) +
            (1 / 24) * c4 * dt^Val(4)
        )

        @test (@ballocated taylor_horner_integral($dt, $cs, $c0)) == 0

        d0 = dimensionless(2.3)

        @test_throws MethodError taylor_horner_integral(dt, cs, d0)
        @test_throws AssertionError taylor_horner(d0, cs)
    end

    @testset "linear algebra" begin
        @testset "dot" begin
            x = [1.1, 0.85, -1.2]
            y = [1.3, -0.87, 1.3]

            sx = distance.(x)
            sy = distance.(y)

            x_y = dot(x, y)
            sx_sy = dot(sx, sy)
            x_sy = dot(x, sy)
            sx_y = dot(sx, y)
            sxTsy = transpose(sx) * sy

            @test x_y ≈ sx_sy.x
            @test x_y ≈ sxTsy.x
            @test x_sy ≈ sx_y
            @test x_sy.x ≈ x_y
            @test udim(sx_sy) == udim(sx[1]) + udim(sy[1])

            @test (@ballocated dot($sx, $sy)) == 0

            x = (1.1, 0.85, -1.2)
            y = (1.3, -0.87, 1.3)

            sx = distance.(x)
            sy = distance.(y)

            x_y = dot(x, y)
            sx_sy = dot(sx, sy)

            @test x_y ≈ sx_sy.x
            @test udim(sx_sy) == udim(sx[1]) + udim(sy[1])

            @constinferred dot(sx, sy)

            @test @ballocated(dot($sx, $sy)) == 0

            @test transpose(x[1]) == x[1]
        end
    end

    @testset "derivatives" begin
        d1 = dimensionless(1.2)
        d2 = dimensionless(2.1)
        d3 = dimensionless(0.5)
        t1 = time(1.3)
        t2 = time(2.4)
        x = 0.9

        @test gradient(a -> value(oneunit(a)), d1) == (dimensionless(0.0),)
        @test @ballocated(
            gradient($(a -> value(oneunit(a))), $d1) == ($dimensionless(0.0),)
        ) == 0

        @test gradient(a -> value(zero(a)), d1) == (dimensionless(0.0),)
        @test @ballocated(
            gradient($(a -> value(zero(a))), $d1) == ($dimensionless(0.0),)
        ) == 0

        @test gradient((t, x) -> value(oftype(t, x)), t1, x) == (zero(1 / t1), 1)
        @test @ballocated(gradient((t, x) -> value(oftype(t, x)), $t1, $x)) == 0

        @test gradient(d -> value(d * d - 2 * d + 1), d1) == (2 * d1 - 2,)
        @test @ballocated(gradient(d -> value(d * d - 2 * d + 1), $d1)) == 0

        @test gradient(d -> value(-d + exp(d)), d1) == (-oneunit(d1) + exp(d1),)
        @test @ballocated(gradient(d -> value(-d + exp(d)), $d1)) == 0

        @test gradient(d -> value((d + 1) / d), d1) == gradient(d -> value(1 + 1 / d), d1)
        @test gradient(d -> value((d + 1) / d), d1)[1] ≈ -1 / d1^2
        @test @ballocated(gradient(d -> value((d + 1) / d), $d1)) == 0
        @test @ballocated(gradient(d -> value(1 + 1 / d), $d1)) == 0

        @test gradient((a, b) -> value(a^b), d1, d2) ==
              gradient((a, b) -> value(exp(b * log(a))), d1, d2)
        @test gradient(a -> value(exp10(a)), d1)[1] ≈ gradient(a -> value(exp(a * log(10))), d1)[1]
        @test gradient(a -> value(exp2(a)), d1)[1] ≈ gradient(a -> value(exp(a * log(2))), d1)[1]
        @test @ballocated(gradient((a, b) -> value(a^b), $d1, $d2)) == 0
        @test @ballocated(gradient((a, b) -> value(exp(b * log(a))), $d1, $d2)) == 0

        @test gradient(a -> value(a^2.0), d1) == gradient(a -> value(a * a), d1)
        @test @ballocated(gradient(a -> value(a^2.0), $d1)) == 0

        @test gradient(a -> value(a^2), d1) == gradient(a -> value(a * a), d1)
        @test_broken @ballocated(gradient(a -> value(a^2), $d1)) == 0

        @test gradient(a -> value(a^unsigned(2)), d1) == gradient(a -> value(a * a), d1)
        @test_broken @ballocated(gradient(a -> value(a^unsigned(2)), $d1)) == 0

        @test gradient(a -> value(a^(3 // 1)), d1) == gradient(a -> value(a * a * a), d1)
        @test @ballocated(gradient(a -> value(a^(3 // 1)), $d1)) == 0

        @test gradient(a -> value(sqrt(a * a * a)), d1) ==
              gradient(a -> value(a^(3 // 2)), d1)
        @test @ballocated(gradient(a -> value(sqrt(a * a * a)), $d1)) == 0

        @test gradient(a -> value(cbrt(sqrt(a))), d1) ==
              gradient(a -> value(root(a, 6)), d1)
        @test @ballocated(gradient(a -> value(cbrt(sqrt(a))), $d1)) == 0
        @test @ballocated(gradient(a -> value(root(a, 6)), $d1)) == 0

        @test gradient(a -> value(root(a, 2 // 3)), d1) ==
              gradient(a -> value(sqrt(a * a * a)), d1)
        @test @ballocated(gradient(a -> value(root(a, 2 // 3)), $d1)) == 0

        @test gradient(a -> value(log2(a)), d1) == gradient(a -> value(log(a) / log(2)), d1)
        @test @ballocated(gradient(a -> value(log2(a)), $d1)) == 0

        @test gradient(a -> value(log10(2 * a)), d1) ==
              gradient(a -> value(log(a) / log(10)), d1)
        @test @ballocated(gradient(a -> value(log10(2 * a)), $d1)) == 0

        @test gradient(a -> value(tan(a)), d1)[1] ≈
              gradient(a -> value(sin(a) / cos(a)), d1)[1]
        @test @ballocated(gradient(a -> value(tan(a)), $d1)) == 0
        @test @ballocated(gradient(a -> value(sin(a) / cos(a)), $d1)) == 0

        @test gradient(a -> value(sec(a)), d1)[1] ≈ gradient(a -> value(1 / cos(a)), d1)[1]
        @test @ballocated(gradient(a -> value(sec(a)), $d1)) == 0

        @test gradient(a -> value(csc(a)), d1)[1] ≈ gradient(a -> value(1 / sin(a)), d1)[1]
        @test @ballocated(gradient(a -> value(csc(a)), $d1)) == 0

        @test gradient(a -> value(cot(a)), d1)[1] ≈
              gradient(a -> value(cos(a) / sin(a)), d1)[1]
        @test @ballocated(gradient(a -> value(cot(a)), $d1)) == 0

        @test gradient(a -> value(asin(a)), d3)[1] == -gradient(a -> value(acos(a)), d3)[1]
        @test @ballocated(gradient(a -> value(asin(a)), $d3)) == 0
        @test @ballocated(gradient(a -> value(acos(a)), $d3)) == 0

        @test gradient(a -> value(atan(sin(a) / cos(a))), d1)[1] == oneunit(d1)
        @test gradient(a -> value(atan(sin(a), cos(a))), d1)[1] == oneunit(d1)
        @test @ballocated(gradient(a -> value(atan(sin(a) / cos(a))), $d1)) == 0
        @test @ballocated(gradient(a -> value(atan(sin(a), cos(a))), $d1)) == 0

        @test gradient(a -> value(asin(sin(a))), d3)[1].x ≈ 1
        @test gradient(a -> value(acos(cos(a))), d3)[1].x ≈ 1
        @test gradient(a -> value(atan(tan(a))), d3)[1].x ≈ 1
        @test gradient(a -> value(acsc(csc(a))), d3)[1].x ≈ 1
        @test gradient(a -> value(asec(sec(a))), d3)[1].x ≈ 1
        @test gradient(a -> value(acot(cot(a))), d3)[1].x ≈ 1

        a = dimensionless(1.1)
        (s, c), _back = pullback(sincos, a)
        @test _back((1, 0))[1] == c
        @test _back((0, 1))[1] == -s
        @test @ballocated(gradient(a -> value(sum(sincos(a))), $a)) == 0

        function func1(a, b, c, t)
            qt = time(t)
            qc = dimensionless(c)
            qb = frequency(b)
            qa = GQ{-2}(a)
            return value(qa * qt * qt + qb * qt + qc)
        end

        function func1_grad_anl(a, b, c, t)
            return t * t, t, 1.0, 2 * a * t + b
        end

        @test collect(func1_grad_anl(1.0, 2.0, 3.0, -2.0)) ≈
              collect(gradient(func1, 1.0, 2.0, 3.0, -2.0))

        function func2(a, w, t)
            qa = dimensionless(a)
            qw = frequency(w)
            qt = time(t)
            return [value(qa * sin(qw * qt)), value(qa * cos(qw * qt))]
        end

        function func2_jac_anl(a, w, t)
            return (
                [sin(w * t), cos(w * t)],
                [a * t * cos(w * t), -a * t * sin(w * t)],
                [a * w * cos(w * t), -a * w * sin(w * t)],
            )
        end

        jac1 = jacobian(func2, Float128(1.2), Float128(0.5), Float128(2.3))
        jac2 = func2_jac_anl(Float128(1.2), Float128(0.5), Float128(2.3))
        @test all([j1 ≈ j2 for (j1, j2) in zip(jac1, jac2)])
    end
end
