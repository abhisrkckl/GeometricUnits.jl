using GeometricUnits
using LinearAlgebra
using Quadmath
using StaticArrays
using Test
using Zygote

@testset verbose = true begin
    t1 = time(2.9)
    t2 = quantity_like(t1, 1.2)

    f1 = frequency(0.2)

    d1 = dimensionless(6.0)
    a1 = dimensionless(0.2)

    @testset "unit and zero" begin
        @test value(unit(t1)) == 1 && udim(unit(t1)) == udim(t1)
        @test value(zero(t1)) == 0 && udim(zero(t1)) == udim(t1)

        @test value(unit(f1)) == 1 && udim(unit(f1)) == udim(f1)
        @test value(zero(f1)) == 0 && udim(zero(f1)) == udim(f1)
    end

    @testset "commonly used quantities" begin
        @test dimensionless(2.5) == speed(2.5)
        @test time(2.5) == distance(2.5)
        @test frequency(2.5) == acceleration(2.5)
    end

    @testset "comparison" begin
        @test t1 == t1
        @test !(t1 != t1)
        @test_throws DomainError t1 == f1

        @test t1 ≈ t1
        @test_throws DomainError t1 ≈ f1

        @test 0 == zero(d1)
        @test 0 ≈ zero(d1)

        @test (t1 > t2) || (t1 <= t2)
        @test (t1 < t2) || (t1 >= t2)

        @test (d1 > 0) || (d1 <= 0)
        @test (d1 < 0) || (d1 >= 0)
        @test (0 > d1) || (0 <= d1)
        @test (0 < d1) || (0 >= d1)

        @test_throws DomainError d1 > t1
        @test_throws DomainError d1 >= t1
        @test_throws DomainError d1 <= t1
        @test_throws DomainError d1 < t1

        @test_throws DomainError 0 > t1
        @test_throws DomainError 0 >= t1
        @test_throws DomainError 0 <= t1
        @test_throws DomainError 0 < t1

        @test isfinite(t1)
        @test !isnan(t1)
        @test !isinf(t1)
    end

    @testset "addition and subtraction" begin
        @test (t1 + t2).x == t1.x + t2.x
        @test_throws DomainError t1 + f1

        @test (t1 - t2).x == t1.x - t2.x
        @test_throws DomainError t1 - f1

        @test (t1 - t2).x == (t1 + (-t2)).x
        @test (t1 + t2).x == (t1 - (-t2)).x

        @test t1 + t2 == t2 + t1
        @test t1 - t2 == -(t2 - t1)
        @test t1 - t1 == zero(t1)

        @test 1 + d1 == d1 + 1
        @test 1 - d1 == -(d1 - 1)

        @test_throws DomainError 1 + t1
        @test_throws DomainError 1 - t1
        @test_throws DomainError t1 + 1
        @test_throws DomainError t1 - 1
    end

    @testset "multiplication and division" begin
        @test (t1 * f1).x == t1.x * f1.x
        @test (t1 * f1).d == t1.d + f1.d

        @test t1 * d1 == d1 * t1
        @test t1 * d1 == d1.x * t1
        @test t1 * d1 == t1 * d1.x

        @test t1 + t1 == 2 * t1
        @test t1 + 3.2 * t1 == 4.2 * t1

        @test (t1 / f1).x == t1.x / f1.x
        @test (t1 / f1).d == t1.d - f1.d

        @test t1 / f1 ≈ t1 * (1 / f1)
        @test t1 / t1 == 1

        @test t1 + 1 / f1 ≈ (t1 * f1 + 1) / f1

        @test 1 / (1 / d1) ≈ d1
        @test 2 * (d1 / 2) ≈ d1
    end

    @testset "power and root" begin
        @test t1 * t1 == t1^2
        @test 1 / t1 == t1^(-1)

        @test sqrt(t1^2) ≈ t1
        @test cbrt(t1^3) ≈ t1
        @test root(t1^5, 5) ≈ t1

        @test sqrt(d1) ≈ d1^0.5
        @test sqrt(d1) ≈ d1^(1 // 2)
        @test d1^2 ≈ root(d1, 1 // 2)

        @test 1.0^d1 ≈ 1
        @test d1^0 ≈ 1

        @test d1^d1 == d1.x^d1.x

        @test_throws DomainError t1^(1 // 2)
        @test_throws DomainError t1^6.3
        @test_throws DomainError 0.5^t1
        @test_throws DomainError sqrt(t1)
        @test_throws DomainError cbrt(t1)
        @test_throws DomainError root(t1, 4)
    end

    @testset "exp and log" begin
        @test exp(log(d1)) ≈ d1
        @test 10^(log10(d1)) ≈ d1
        @test 2^(log2(d1)) ≈ d1
    end

    @testset "trigonometric functions" begin
        @test sin(d1) ≈ 1 / csc(d1)
        @test cos(d1) ≈ 1 / sec(d1)
        @test tan(d1) ≈ 1 / cot(d1)

        @test asin(sin(a1)) ≈ a1
        @test acos(cos(a1)) ≈ a1
        @test atan(tan(a1)) ≈ a1
        @test acsc(csc(a1)) ≈ a1
        @test asec(sec(a1)) ≈ a1
        @test acot(cot(a1)) ≈ a1

        @test atan(t1, t2) ≈ atan(t1 / t2)
    end

    @testset "taylor-horner" begin
        t0 = time(0.0)
        c0 = distance(1.5)
        c1 = speed(1.0)
        c2 = acceleration(1.2)
        c3 = acceleration(0.5) / time(2.0)
        c4 = acceleration(1.2) / time(2.0)^2
        cs = [c1, c2, c3, c4]
        th = TaylorSeries(t0, c0, cs)

        t1 = time(1.0)

        @test th(t1) ≈
              c0 +
              c1 * t1 +
              (1 / 2) * c2 * t1^2 +
              (1 / 6) * c3 * t1^3 +
              (1 / 24) * c4 * t1^4

        d0 = dimensionless(2.3)

        @test_throws DomainError th(d0)

        @test_throws DomainError TaylorSeries(t0, d0, cs)

        c5 = acceleration(1.2) / time(2.0)^2
        cs = [c1, c2, c3, c4, c5]
        @test_throws DomainError TaylorSeries(t0, c0, cs)
    end

    @testset "linear algebra" begin
        @testset "dot" begin
            x = [1.1, 0.85, -1.2]
            y = [1.3, -0.87, 1.3]
            M = [[1.0, 2.0, 3] [4, 4.3, -2] [3.0, 2.1, -3.0]]

            sx = distance.(x)
            sy = distance.(y)

            x_y = dot(x, y)
            sx_sy = dot(sx, sy)

            @test x_y ≈ sx_sy.x
            @test sx_sy.d == sx[1].d + sy[1].d

            sax = distance.(x)
            say = distance.(y)
            
            sax_say = dot(sax, say)

            @test x_y ≈ sax_say.x
            @test sax_say.d == sax[1].d + say[1].d
        end
    end

    @testset "derivatives" begin
        function func1(a, b, c, t)
            qt = time(t)
            qc = dimensionless(c)
            qb = frequency(b)
            qa = GQ(a, -2)
            return value(qa * qt * qt + qb * qt + qc)
        end

        function func1_grad_anl(a, b, c, t)
            return t * t, t, 1.0, 2 * a * t + b
        end

        @test [func1_grad_anl(1.0, 2.0, 3.0, -2.0)...] ≈ [gradient(func1, 1.0, 2.0, 3.0, -2.0)...]

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
