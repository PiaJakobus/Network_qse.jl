@testset "prefactor" begin
    a = Network_qse.extract_partition_function()
    find_el = filter(i -> (a[i].name == "He4"), 1:size(a,1))[1]
    ni56 = filter(i -> (a[i].name == "Ni56"), 1:size(a,1))
    he3 = Network_qse.AtomicProperties(0,0, 0.5, 2.39e-5, o -> 1.0)
    o16 = Network_qse.AtomicProperties(8,16, 0.0, -4.737, o -> 1.0)
    @test Network_qse.prefactor(a[find_el])(2e9, 20235896) / 44161 * 100 - 100 < 1
    @test Network_qse.prefactor(he3)(1e9, 1e7) == 0.0
    @test (100 * (1 - 28555018 / Network_qse.prefactor(o16)(2e9,1e6))) < 1
end


@testset "initial_guess" begin
    pf = Network_qse.extract_partition_function()
    @test Network_qse.initial_guess(pf)[1] / -9.17 < 1.01
end


@testset "df_nse_condition" begin
    a = Network_qse.extract_partition_function()
    x = Network_qse.initial_guess(a)
    g(x) = Network_qse.nse_condition(x, 9e9, 1e7, 0.49, a)
    @test Network_qse.df_nse_condition(x, 9e9, 1e7, a) ≈ Network_qse.ForwardDiff.jacobian(g, x)
end


@testset "df_qse_condition" begin
    a = Network_qse.extract_partition_function()
    x = Network_qse.qse_initial_guess(a)
    @test Network_qse.df_qse_condition(x, 9e9, 1e7, a) ≈ Network_qse.ForwardDiff.jacobian(x -> Network_qse.qse_condition(x, 9e9, 1e7, 0.5, 1.0, a), x)
end

@testset "qse_condition" begin
    a = Network_qse.extract_partition_function()
    tmp = [-9.162532665782894, -9.162532665782894, -493.08391464192107]
    @test Network_qse.qse_condition(tmp, 7.25e9, 1e7, 0.5, 0.9, a)[1] ≈ Network_qse.nse_condition(tmp[1:2], 7.25e9, 1e7, 0.5, a)[1]
    @test Network_qse.qse_condition(tmp, 7.25e9, 1e7, 0.5, 0.9, a)[2] ≈ Network_qse.nse_condition(tmp[1:2], 7.25e9, 1e7, 0.5, a)[2]
end
