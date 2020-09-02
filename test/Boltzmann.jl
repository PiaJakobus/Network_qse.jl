@testset "prefactor" begin
    #TODO: more reliable tests for chargeNeut, massCon
    #TODO: why is type declaration for args::AtomicProperties not working?
    he3 = Network_qse.AtomicProperties(0,0, 0.5, 2.39e-5, o -> 1.0)
    @test Network_qse.prefactor(he3)(1e9, 1e7) == 0.0
    #@test Network_qse.chargeNeut([0.0000001,-0.000007], 1e9, 1e7, he3) < 10.0
    #@test Network_qse.massCon([0.0000001,-0.000007], 1e9, 1e7, he3) < 10.0
end
