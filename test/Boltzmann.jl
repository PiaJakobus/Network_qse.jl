@testset "prefactor" begin
    #TODO: why is type declaration for args::AtomicProperties not working?
    he3 = Network_qse.AtomicProperties(0,0, 0.5, 2.39e-5, o -> 1.0)
    o16 = Network_qse.AtomicProperties(8,16, 0.0, -4.737, o -> 1.0)
    @test Network_qse.prefactor(he3)(1e9, 1e7) == 0.0
    @test (100 * (1 - 28555018 / Network_qse.prefactor(o16)(2e9,1e6))) < 1
    #@test Network_qse.chargeNeut([0.0000001,-0.000007], 1e9, 1e7, he3) < 10.0
end


@testset "initial_guess" begin
    pf = Network_qse.extract_partition_function()
    ind = filter(i -> (pf[i].name == "Fe56"), 1:size(pf,1))
    @test Network_qse.initial_guess(2e9, 1e6, pf[ind][1]) / -9.17 < 1.01
end
