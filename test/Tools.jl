@testset "logsumexp" begin
    @test Network_qse.logsumexp([1.0,2.0,3.0]) ≈ log(sum(exp(1.0) + exp(2.0) + exp(3.0)))
end

@testset "MultiNewtonRaphson" begin
    a = Network_qse.extract_partition_function()
    @test Network_qse.MultiNewtonRaphson([-9.2,-9.3], 9e9, 1e7, 0.5, a) ≠ NaN
end

@testset "inv_3x3" begin
    @test Network_qse.inv_3x3(Matrix(1.0I, 3, 3)) == Matrix(1.0I, 3, 3)
end

#TODO: add test
@testset "QSE_MultiNewtonRaphson" begin
    1 == 1
end
