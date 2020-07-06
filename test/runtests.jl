using Test
using Network_qse


@testset "network_nse" begin
    @test isa(initial_partition_function()[1], Array{Float64,2})
end

@testset "tools" begin
    @test findnearest([1:10;], 10) == 10:10
    @test linear_interpolation([1.,2.], [1.,2.], 10.) == 10.
end

@testset "io" begin
    @test isa(read_part_frdm(), Tuple{Array{Float64,2}, Array{Array{Float64,2},1}})
    @test isa(read_species(), Tuple{Array{Float64,2},Int64})
    @test read_mass_frdm()[66,:] == [7.0, 13.0, 2.0, 5.345]
    @test isa(extract_partition_function()[2:5], NTuple{4,Array{Float64,1}})
end
