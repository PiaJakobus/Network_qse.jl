@testset "prefactor" begin
    #TODO: more reliable tests for chargeNeut, massCon
    #TODO: why is type declaration for args::AtomicProperties not working?
    he3 = Network_qse.AtomicProperties(0,0, 0.5, 2.39e-5, o -> 1.0)
    @test Network_qse.prefactor(he3)(1e9, 1e7) == 0.0
    #@test Network_qse.chargeNeut([0.0000001,-0.000007], 1e9, 1e7, he3) < 10.0
    #@test Network_qse.massCon([0.0000001,-0.000007], 1e9, 1e7, he3) < 10.0
end


@testset "initial_guess" begin
    pf = Network_qse.extract_partition_function()
    ind = filter(i -> (pf[i].name == "Fe56"), 1:size(pf,1))
    @test Network_qse.initial_guess(2e9, 1e6, pf[ind][1]) / -9.17 < 1.01  
end
