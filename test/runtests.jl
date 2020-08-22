using Test
using Network_qse


@testset "network_nse" begin
    ω,A,Z,s,m = extract_partition_function()
    @test isa(initial_partition_function(ω,A,Z,s,m), Array{BigFloat,2})
    @test logsumexp([1,2,3]) - log(sum(exp.([1,2,3]))) < 1e-10
    @test eos((1,2,3))[1][1] ==  1.2067926406393289e6
    @test sum(mass_fraction([1.5e-5,1.1e-5], 3e9, 1e9,A,Z,m)) - exp(logsumexp(log_mass_fraction([1.5e-5,1.1e-5], 3e9, 1e9,A,Z,m))) < 1e-5
    @test find_nucl(56,26,A,Z) = 661
end

@testset "tools" begin
    @test findnearest([1:10;], 10) == 10:10
    @test linear_interpolation([1.,2.], [1.,2.], 10.) == 10.
end

@testset "io" begin
    @test isa(read_part_frdm()[1], Array{BigFloat,2})
    @test isa(read_species(), Tuple{Array{BigFloat,2},Int64})
    @test isa(read_mass_frdm()[4,:],Array{BigFloat,1})
    @test isa(extract_partition_function()[2:5], NTuple{4})
end

@testset "DataTypes" begin
    include("DataTypes.jl")
end

@testset "Boltzmann" begin
    include("Boltzmann.jl")
end
