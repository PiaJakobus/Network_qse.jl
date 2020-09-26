@testset "Tools" begin
    a = Network_qse.extract_partition_function()
    @test Network_qse.logsumexp([1.0,2.0,3.0]) ≈ log(sum(exp(1.0) + exp(2.0) + exp(3.0)))
    @test Network_qse.MultiNewtonRaphson([-9.2,-9.3], 9e9, 1e7, 0.5, a) ≠ NaN
end
