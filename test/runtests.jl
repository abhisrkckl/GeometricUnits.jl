using GeometricUnits
using Test

@testset verbose = true begin
    t1 = time(2.9)
    t2 = quantity(t1, 1.2)

    f1 = frequency(0.2)

    d1 = dimensionless(6.0)

    @testset "unit and zero" begin
        @test unit(t1).x == 1 && unit(t1).d == t1.d
        @test zero(t1).x == 0 && zero(t1).d == t1.d

        @test unit(f1).x == 1 && unit(f1).d == f1.d
        @test zero(f1).x == 0 && zero(f1).d == f1.d
    end

    @testset "comparison" begin
        @test t1 == t1
        @test_throws DomainError t1 == f1

        @test t1 ≈ t1
        @test_throws DomainError t1 ≈ f1
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

        @test t1 + 1 / f1 ≈ (t1 * f1 + 1)
    end
end
