@testset "prefactor" begin
    k = AtomicProperties(3,3,1.0,10.0)
    pf1 = PartitionFunction(k,1.0,1e9)
    pf2 = PartitionFunction(k,1.0,2e9)
    l   = [pf1,pf2]
    @test prefactor(pf1) ≈ 4.637634538831878e11
    @test all(prefactor.(l) .≈ [4.637634538831878e11,1.3117211324291873e12])
end
