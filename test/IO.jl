@testset "IO" begin
    g    = Network_qse.extract_partition_function()
    fe56 = filter(i -> g[i].name == "Fe56", 1:size(g,1))[1]
    p   = filter(i -> g[i].name == "H1", 1:size(g,1))[1]
    H2   = filter(i -> g[i].name == "H2", 1:size(g,1))[1]
    n    = filter(i -> g[i].name == "1n", 1:size(g,1))[1]
    th   = Network_qse.ThermoProperties(1e9, 1e7, 0.5, -10.0)
    @test isa(Network_qse.read_part_frdm()[1], Array{Float64,2})
    @test isa(Network_qse.read_species(), Tuple{Array{Float64,2},Int64})
    @test isa(Network_qse.read_mass_frdm()[4,:],Array{Float64,1})
    @test isa(g, Array{Network_qse.AtomicProperties,1})
    @test g[n].ω(1e9) == 2.0
    @test g[n].Z == 0
    @test g[n].A == 1
    @test g[n].M / (Network_qse.m_n*Network_qse.ergmev*Network_qse.c^2) < 1.02
    @test g[n].s == 0.5
    @test g[p].Z == 1
    @test g[p].A == 1
    @test g[p].M / (Network_qse.m_n*Network_qse.ergmev*Network_qse.c^2) < 1.02
    @test g[H2].A == 2
    @test g[H2].Z == 1
    @test g[fe56].name == "Fe56"
    @test g[fe56].M / (9.29e-23*Network_qse.ergmev*Network_qse.c^2) < 1.02
    @test g[fe56].Δ / (9.708e-5*Network_qse.ergmev*Network_qse.c^2) < 1.02
    @test (g[fe56].Eb/56) / 8.8 < 1.02
    @test th.T == 1e9 && th.rho == 1e7 && th.y == 0.5 && th.x_qse == 0.9999999999
end
