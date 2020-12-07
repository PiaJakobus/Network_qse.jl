@testset "logsumexp" begin
    @test Network_qse.logsumexp([1.0,2.0,3.0]) ≈ log(sum(exp(1.0) + exp(2.0) + exp(3.0)))
end

@testset "MultiNewtonRaphson" begin
    a = Network_qse.extract_partition_function()
    x = Network_qse.initial_guess(a)
    th = Network_qse.ThermoProperties(9e9, 1e7, 0.49, -5.0)
    sol = Network_qse.MultiNewtonRaphson(x, th, a)
    @test round(sum(Network_qse.x_i(sol, th, a)), digits = 1) ≈ 1.0
    @test Network_qse.MultiNewtonRaphson([-9.2,-9.3],th, a) ≠ NaN
end

@testset "inv_3x3" begin
    @test Network_qse.inv_3x3(Matrix(1.0I, 3, 3)) == Matrix(1.0I, 3, 3)
end

@testset "find_el" begin
    a = Network_qse.extract_partition_function()
    @test a[Network_qse.find_el("1n", a)].Eb == 0.0
    @test a[Network_qse.find_el("H1", a)].Eb == 0.0
    @test a[Network_qse.find_el("Fe56", a)].Eb ≈ -492.245
end

#TODO: add test
@testset "QSE_MultiNewtonRaphson" begin
    1 == 1
end
