@testset "prefactor" begin
    #TODO: why is type declaration for args::AtomicProperties not working?
    a = Network_qse.extract_partition_function()
    find_el = filter(i -> (a[i].name == "He4"), 1:size(a,1))[1]
    ni56 = filter(i -> (a[i].name == "Ni56"), 1:size(a,1))
    he3 = Network_qse.AtomicProperties(0,0, 0.5, 2.39e-5, o -> 1.0)
    o16 = Network_qse.AtomicProperties(8,16, 0.0, -4.737, o -> 1.0)
    @test Network_qse.prefactor(a[find_el])(2e9, 20235896) / 44161 * 100 - 100 < 1
    @test Network_qse.prefactor(he3)(1e9, 1e7) == 0.0
    @test (100 * (1 - 28555018 / Network_qse.prefactor(o16)(2e9,1e6))) < 1
    #@test Network_qse.chargeNeut([0.0000001,-0.000007], 1e9, 1e7, he3) < 10.0
end


@testset "initial_guess" begin
    pf = Network_qse.extract_partition_function()
    ind = filter(i -> (pf[i].name == "Fe56"), 1:size(pf,1))
    @test Network_qse.initial_guess(pf[ind][1])[1] / -9.17 < 1.01
end


@testset "df_nse_condition" begin
    a = Network_qse.extract_partition_function()
    x = [-9.1, -9.1]
    sol = Network_qse.MultiNewtonRaphson([-9.2,-9.3], 9e9, 1e9, 0.49, a)
    g(x) = Network_qse.nse_condition(x, 3e9, 1e7, 0.49, a)
    @test all(Network_qse.df_nse_condition(x, 3e9, 1e7, 0.49, a) .≈ Network_qse.ForwardDiff.jacobian(g, x))
    @test sum(Network_qse.x_i(sol, 9e9, 1e9, a)) ≈ 1.0
end
