@testset "AtomicProperties" begin
    t1 = AtomicProperties(2,3,1.0,10.0)
    @test t1.name == "He3"
    # TODO: more tests!
end

@testset "PartitionFunction" begin
    k = AtomicProperties(3,3,1.0,10.0)
    @test k.T == 1.0e9
    #TODO: more tests!
end
