@testset "logsumexp" begin
    @test Network_qse.logsumexp([1.0,2.0,3.0]) ≈ log(sum(exp(1.0) + exp(2.0) + exp(3.0)))
end

@testset "MultiNewtonRaphson" begin
    a = Network_qse.extract_partition_function()
    x = Network_qse.qse_initial_guess(a)
    x1 = Network_qse.initial_guess(a)
    th = Network_qse.ThermoProperties(5e9, 1e7, 0.49, -2.0)
    th1 = Network_qse.ThermoProperties(5e9, 1e7, 0.49, -2.0)
    sp1 =  Network_qse.StepParameter(-50,50,30)
    sp =  Network_qse.StepParameter(-10,10,50)

    ff = Network_qse.Func(3, x -> Network_qse.qse_condition(x, th, a), x -> Network_qse.df_qse_condition(x,th,a), false)
    ff1 = Network_qse.Func(2, x -> Network_qse.nse_condition(x, th, a), x -> Network_qse.df_nse_condition(x,th,a), false)
    sol = Network_qse.MultiNewtonRaphson(x, ff, th, a, sp)
    sol1 = Network_qse.MultiNewtonRaphson(x1, ff1, th, a, sp1)
    println(sol,sol1)
    @test sum(Network_qse.x_i_QSE(sol, th, a)[2]) ≈ th.x_qse
    @test sum(Network_qse.x_i(sol1, th1, a)) ≈ 1.0
end

@testset "inv_3x3" begin
    @test Network_qse.inv_3x3(Matrix(1.0I, 3, 3)) == Matrix(1.0I, 3, 3)
end

@testset "inv_2x2" begin
    @test 1 == 1
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
