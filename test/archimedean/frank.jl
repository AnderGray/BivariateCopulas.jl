@testset "Frank" begin
    @testset "constructor" begin
        @suppress @test isa(Frank(0), Independence)
        @test_throws AssertionError Frank(-2.0)

        @test_logs (:warn, "Frank returns an independence copula for ϑ == 0.0") Frank(
            0.0
        )
        @test_logs (:warn, "Frank returns an M copula for ϑ == Inf") Frank(Inf)
    end
    @testset "τ" begin
        @test τ(Frank(1.0)) == 0.11001853644899295
        @test τ(Frank(10.0)) == 0.6657773862719785
        @test τ(Frank(20.0)) == 0.81644934023564
    end

end
