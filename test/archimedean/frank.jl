n = 10^5
θ = [1.0, 10.0, 20.0]

u = collect(range(0.0, 1.0, length=10))
v = u

@testset "Frank" begin
    @testset "constructor" begin
        @test_throws AssertionError Frank(Inf)
        @test isa(Frank(0), Independence)

        @test_logs (:warn, "Frank returns an independence copula for ϑ == 0.0") Frank(
            0.0
        )
    end

    @testset "generators" begin
        u = [0:0.1:1;]
        c = Frank(10.0)
        @test BivariateCopulas.φ⁻¹.(BivariateCopulas.φ.(u, c), c) ≈ u
    end

    @testset "τ" begin
        @test τ(Frank(1.0)) == 0.11001853644899295
        @test τ(Frank(10.0)) == 0.6657773862719785
        @test τ(Frank(20.0)) == 0.81644934023564
    end

    @testset "sample" begin
        c = Frank(10.0)
        u = BivariateCopulas.sample(c, n)
        @test corkendall(u) ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01
    end

    @testset "rosenblatt" begin
        c = Frank(20.0)
        u = BivariateCopulas.sample(c, n)
        @test corkendall(rosenblatt(u, c)) ≈ [1.0 0.0; 0.0 1.0] atol = 0.01
    end

    @testset "inverse_rosenblatt" begin
        c = Frank(1.0)
        u = BivariateCopulas.sample(c, n)

        @test inverse_rosenblatt(rosenblatt(u, c), c) ≈ u
    end

    @testset "sklar's theorem" begin
        c = Frank(1.0)
        X = Normal()
        Y = Normal()

        @test isa(c(X, Y), Joint)
    end

    @testset "cdf" begin
        x = [0:0.25:1;]
        y = x

        @test cdf.(Frank(10.0), x, y) ≈
              [0.0, 0.18490043593650232, 0.4313568167929175, 0.6849004359365126, 0.9999999999999871]
        @test cdf.(Frank(1.0), x, y) ≈
              [0.0, 0.08056458733770672, 0.28092980362016146, 0.580564587337707, 1.0]
    end

    @testset "density" begin
        x = [0:0.25:1;]
        y = x

        @test isnan(BivariateCopulas.density.(Frank(10.0), x[1], y[1]))
        @test BivariateCopulas.density.(Frank(10.0), x[2:end], y[2:end]) ≈
              [2.720019944137354, 2.533918274531531, 2.7200199441379125, 10.000454019907506]

        @test isnan(BivariateCopulas.density.(Frank(1.0), x[1], y[1]))
        @test BivariateCopulas.density.(Frank(1.0), x[2:end], y[2:end]) ≈
              [1.127276244650801, 1.0207470412683992, 1.1272762446508016, 1.5819767068693265]
    end

    @testset "grounded" begin
        for ϑ in θ
            c = Frank(ϑ)
            @test iszero(cdf.(c, u, 0.0))
            @test iszero(cdf.(c, 0.0, v))
        end
    end

    @testset "uniform margins" begin
        for ϑ in θ
            c = Frank(ϑ)
            @test cdf.(c, u, 1.0) ≈ u
            @test cdf.(c, 1.0, v) ≈ v
        end
    end

    @testset "2-increasing" begin
        for ϑ in θ
            c = Frank(ϑ)
            for i in 1:10
                for j in i:10
                    @test cdf(c, u[j], v[j]) - cdf(c, u[j], v[i]) - cdf(c, u[i], v[j]) + cdf(c, u[i], v[i]) >= 0.0
                end
            end
        end
    end
end
