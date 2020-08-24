@testset "Tools" begin
    @test Network_qse.logsumexp([1.0,2.0,3.0]) ≈ log(sum(exp(1.0) + exp(2.0) + exp(3.0)))
end
