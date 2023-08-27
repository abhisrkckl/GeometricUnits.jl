using GeometricUnits
using Test

@testset verbose = true begin
    t1 = time(2.9)
    t2 = quantity_like(t1, 1.2)

    f1 = frequency(0.2)

    d1 = dimensionless(6.0)
    a1 = dimensionless(0.2)

    @testset "unit and zero" begin
        @test unit(t1).x == 1 && unit(t1).d == t1.d
        @test zero(t1).x == 0 && zero(t1).d == t1.d

        @test unit(f1).x == 1 && unit(f1).d == f1.d
        @test zero(f1).x == 0 && zero(f1).d == f1.d
    end

    @testset "commonly used quantities" begin
        @test dimensionless(2.5) == speed(2.5)
        @test time(2.5) == distance(2.5)
        @test frequency(2.5) == acceleration(2.5)
    end

    @testset "comparison" begin
        @test t1 == t1
        @test_throws DomainError t1 == f1

        @test t1 ≈ t1
        @test_throws DomainError t1 ≈ f1

        @test 0 == zero(d1)
        @test 0 ≈ zero(d1)
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
end
