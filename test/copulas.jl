copulas = [Clayton.([-0.75, 2])... Frank.([1, 10])...]
n = 10^5

@testset "Copulas" begin
    @testset "grounded" begin
        u = range(0.0, 1.0, length=10)
        v = u
        for c in copulas
            @test iszero(cdf.(c, u, 0.0))
            @test iszero(cdf.(c, 0.0, v))
        end
    end

    @testset "uniform margins" begin
        u = range(0.0, 1.0, length=10)
        v = u
        for c in copulas
            @test cdf.(c, u, 1.0) ≈ u
            @test cdf.(c, 1.0, v) ≈ v
        end
    end

    @testset "2-increasing" begin
        u = range(0.0, 1.0, length=100)
        v = u
        for c in copulas
            for i in 1:100
                for j in i:100
                    @test cdf(c, u[j], v[j]) - cdf(c, u[j], v[i]) - cdf(c, u[i], v[j]) + cdf(c, u[i], v[i]) >= 0.0
                end
            end
        end
    end

    @testset "sample" begin
        for c in copulas
            u = BivariateCopulas.sample(c, n)
            @test corkendall(u) ≈ [1.0 τ(c); τ(c) 1.0] atol = 0.01
        end
    end

    @testset "sklar's theorem" begin
        X = Normal()
        Y = Normal()

        for c in copulas
            @test isa(c(X, Y), Joint)
        end
    end

    @testset "rosenblatt" begin
        for c in copulas
            u = BivariateCopulas.sample(c, n)
            @test corkendall(rosenblatt(u, c)) ≈ [1.0 0.0; 0.0 1.0] atol = 0.01
        end
    end

    @testset "inverse rosenblatt" begin
        for c in copulas
            u = BivariateCopulas.sample(c, n)

            @test inverse_rosenblatt(rosenblatt(u, c), c) ≈ u
        end
    end

    @testset "density" begin
        for c in copulas

            if isa(c, Clayton) && c.ϑ < 0 # takes too long to converge otherwise
                0.95 <= hcubature(x -> density(c, x[1], x[2]), zeros(2), ones(2); maxevals=10^6)[1] <= 1.05
                continue
            end
            A = hcubature(x -> density(c, x[1], x[2]), zeros(2), ones(2))
            @test A[1] ≈ 1.0 atol = A[2]

            x = range(0, 1; length=10)
            y = range(0, 1; length=10)
            for xᵢ in x, yᵢ in y
                @test density(c, xᵢ, yᵢ) >= 0
            end
        end
    end
end
