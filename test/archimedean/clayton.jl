n = 10^6
θ = [-0.5, 2, 10]

@testset "Clayton" begin
    @testset "constructor" begin
        @test isa(Clayton(0), Independence)
        @test_throws AssertionError Clayton(-2)

        @test_logs (:warn, "Clayton returns a W copula for ϑ < -0.5") Clayton(-0.7)
        @test_logs (:warn, "Clayton returns an independence copula for ϑ == 0.0") Clayton(
            0.0
        )
        @test_logs (:warn, "Clayton returns an M copula for ϑ == Inf") Clayton(Inf)
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
            @test corkendall(u) ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01
        end
    end

    @testset "rosenblatt" begin
        for ϑ in θ
            c = Clayton(ϑ)
            u = BivariateCopulas.sample(c, n)
            @test corkendall(rosenblatt(u, c)) ≈ [1.0 0.0; 0.0 1.0] atol = 0.01
        end
    end

    @testset "inverse_rosenblatt" begin
        for ϑ in θ
            c = Clayton(ϑ)
            u = BivariateCopulas.sample(c, n)

            @test inverse_rosenblatt(rosenblatt(u, c), c) ≈ u
        end
    end

    @testset "sklar's theorem" begin
        c = Clayton(2)
        X = Normal()
        Y = Normal()

        @test isa(c(X, Y), Joint)
    end

    # TODO: Test density
end
