@testset "Clayton" begin
    @testset "constructor" begin
        @suppress @test isa(Clayton(0), Independence)
        @test_throws AssertionError Clayton(-2)

        @test_logs (:warn, "Clayton returns an independence copula for ϑ == 0.0") Clayton(
            0.0
        )
        @test_logs (:warn, "Clayton returns an M copula for ϑ == Inf") Clayton(Inf)
        @test_logs (:warn, "Clayton returns a W copula for ϑ == -1") Clayton(-1)
    end

    @testset "τ" begin
        @test τ(Clayton(-0.75)) == -0.6
        @test τ(Clayton(1)) == 1 / 3
        @test τ(Clayton(2)) == 0.5
    end

end
