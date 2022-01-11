n = 10^6
θ = [-0.5, 2, 10]

function _clayton_density(u, v, ϑ)
    return (ϑ + 1) * (u * v)^(-ϑ - 1) * (u^(-ϑ) + v^(-ϑ) - 1)^(-(2ϑ + 1) / ϑ)
end

@testset "Clayton" begin
    @testset "constructor" begin
        @test isa(Clayton(0), Independence)
        @test_throws AssertionError Clayton(-1)
        @test_throws AssertionError Clayton(-2)
        @test_throws AssertionError Clayton(Inf)
    end

    @testset "τ" begin
        @test τ(Clayton(θ[1])) == -1 / 3
        @test τ(Clayton(θ[2])) == 0.5
        @test τ(Clayton(θ[3])) == 10 / 12
    end

    @testset "sample" begin
        for ϑ in θ
            c = Clayton(ϑ)
            u = BivariateCopulas.sample(c, n)
            @test corkendall(u') ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01
        end
    end

    @testset "rosenblatt" begin
        for ϑ in θ
            c = Clayton(ϑ)
            u = BivariateCopulas.sample(c, n)
            @test corkendall(rosenblatt(u, c)') ≈ [1.0 0.0; 0.0 1.0] atol = 0.01
        end
    end

    @testset "inverse_rosenblatt" begin
        for ϑ in θ
            c = Clayton(ϑ)
            u = BivariateCopulas.sample(c, n)

            @test inverse_rosenblatt(rosenblatt(u, c), c) ≈ u
        end
    end

    @testset "density" begin
        u = range(0.01, 1, 10)
        v = range(0.01, 1, 10)

        for ϑ in θ
            @test density.(Clayton(ϑ), u, v) ≈ _clayton_density.(u, v, ϑ)
        end
    end
end
