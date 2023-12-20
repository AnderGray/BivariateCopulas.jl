copulas = [Clayton.([-0.75, 1])... Frank.([1.0, 10.0])...]

@testset "Archimedean Copulas" begin

    @testset "Generators" begin
        u = range(0, 10, length=10)
        for c in copulas
            @test BivariateCopulas.φ⁻¹.(BivariateCopulas.φ.(u, c), c) ≈ u
            @test BivariateCopulas.φ(0.0, c) ≈ 1.0

            if isa(c, Clayton) && c.ϑ < 0 # not strict for ϑ < 0
                @test BivariateCopulas.φ(floatmax(), c) == Inf
                @test issorted(BivariateCopulas.φ.(u, c))
            else
                @test BivariateCopulas.φ(floatmax(), c) ≈ 0.0 atol = 1e-20
                @test issorted(BivariateCopulas.φ.(u, c); rev=true)
            end
        end
    end

    @testset "Inverse Generators" begin
        for c in copulas
            if isa(c, Clayton) && c.ϑ < 0 # not strict for ϑ < 0
                @test BivariateCopulas.φ⁻¹(0.0, c) == -1.0
            else
                @test BivariateCopulas.φ⁻¹(0.0, c) == Inf
            end
            @test BivariateCopulas.φ⁻¹(1.0, c) == 0.0
        end
    end
end
