n = 10^5
θ = [1, 2, 10]

@testset "Clayton" begin
    @testset "constructor" begin
        @test isa(Clayton(0), Independence)
        @test_throws AssertionError Clayton(-2)

        @test_logs (:warn, "Clayton returns an independence copula for ϑ == 0.0") Clayton(
            0.0
        )
        @test_logs (:warn, "Clayton returns an M copula for ϑ == Inf") Clayton(Inf)
    end

    @testset "generators" begin
        u = [0:0.1:1;]
        for ϑ in θ
            c = Clayton(ϑ)
            @test BivariateCopulas.φ⁻¹.(BivariateCopulas.φ.(u, c), c) ≈ u
        end
    end

    @testset "τ" begin
        @test τ(Clayton(θ[1])) == 1 / 3
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

    @testset "cdf" begin
        x = [0:0.25:1;]
        y = x

        @test cdf.(Clayton(2), x, y) ≈
              [0.0, 0.1796053020267749, 0.37796447300922725, 0.6255432421712244, 1.0]
    end

    @testset "density" begin
        x = [0:0.25:1;]
        y = x

        @test isnan(BivariateCopulas.density(Clayton(2), x[1], y[1]))
        @test BivariateCopulas.density.(Clayton(2), x[2:end], y[2:end]) ≈
              [2.2965556205046926, 1.481003649342278, 1.614508582188617, 3.0]
    end

    @testset "grounded" begin
        u = collect(range(0.0, 1.0, 10))
        v = u
        for ϑ in θ
            c = Clayton(ϑ)
            @test iszero(cdf.(c, u, 0.0))
            @test iszero(cdf.(c, 0.0, v))
        end
    end

    @testset "uniform margins" begin
        u = collect(range(0.0, 1.0, 10))
        v = u
        for ϑ in θ
            c = Clayton(ϑ)
            @test cdf.(c, u, 1.0) ≈ u
            @test cdf.(c, 1.0, v) ≈ v
        end
    end

    @testset "2-increasing" begin
        u = collect(range(0.0, 1.0, 10))
        v = u
        for ϑ in θ
            c = Clayton(ϑ)
            for i in 1:10
                for j in i:10
                    @test cdf(c, u[j], v[j]) - cdf(c, u[j], v[i]) - cdf(c, u[i], v[j]) + cdf(c, u[i], v[i]) >= 0.0
                end
            end
        end
    end
end
