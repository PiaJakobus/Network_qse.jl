@testset "IO" begin
    g = Network_qse.extract_partition_function()
    @test isa(Network_qse.read_part_frdm()[1], Array{Float64,2})
    @test isa(Network_qse.read_species(), Tuple{Array{Float64,2},Int64})
    @test isa(Network_qse.read_mass_frdm()[4,:],Array{Float64,1})
    @test isa(g, Array{Network_qse.AtomicProperties,1})
    @test g[1].ω(1e9) == 2.0
    @test g[1].Z == 0
    @test g[1].A == 1
    @test g[1].m == Network_qse.m_n
    @test g[1].s == 0.5
    @test g[2].Z == 1
    @test g[2].A == 1
    @test g[2].m == Network_qse.m_p
    @test g[3].A == 2
    @test g[3].Z == 1
end
