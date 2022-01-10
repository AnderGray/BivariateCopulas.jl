n = 10^6

@testset "Clayton" begin
    @testset "constructor" begin
        @test isa(Clayton(0), Independence)
        @test_throws AssertionError Clayton(-1)
        @test_throws AssertionError Clayton(-2)
        @test_throws AssertionError Clayton(Inf)
    end

    @testset "τ" begin
        @test τ(Clayton(-1 / 2)) == -1 / 3
        @test τ(Clayton(2)) == 0.5
        @test τ(Clayton(10)) == 10 / 12
    end

    @testset "sample" begin
        c = Clayton(-1 / 2)
        u = BivariateCopulas.sample(c, n)
        @test corkendall(u') ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01

        c = Clayton(2)
        u = BivariateCopulas.sample(c, n)
        @test corkendall(u') ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01

        c = Clayton(10)
        u = BivariateCopulas.sample(c, n)
        @test corkendall(u') ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01
    end

    @testset "rosenblatt" begin
        c = Clayton(-0.5)
        u = BivariateCopulas.sample(c, n)

        @test corkendall(rosenblatt(u, c)') ≈ [1.0 0.0; 0.0 1.0] atol = 0.01

        c = Clayton(2)
        u = BivariateCopulas.sample(c, n)

        @test corkendall(rosenblatt(u, c)') ≈ [1.0 0.0; 0.0 1.0] atol = 0.01

        c = Clayton(10)
        u = BivariateCopulas.sample(c, n)

        @test corkendall(rosenblatt(u, c)') ≈ [1.0 0.0; 0.0 1.0] atol = 0.01
    end
end
