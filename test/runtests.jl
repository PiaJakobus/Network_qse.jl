using Test
using Network_qse
#using DelimitedFiles


@testset "network_nse" begin
    @test Network_qse.testing() == 1
    #@test isa(initial_partition_function(1,1), Array{Float64,1})
end

@testset "io" begin
    @test isa(read_part_frdm(), Tuple{Array{Float64,2},Array{Array{Float64,1},2}})
    #@test isa(read_mass_frdm(), Array{Float64,2})
    @test isa(read_species(), Tuple{Array{Float64,2},Int64})
    @test isa(extract_partition_function(), Tuple{Array{Array{Float64,2},1},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}})
    @test findnearest([1:10;], 10) == 10:10
    @test linear_interpolation([1.,2.], [1.,2.], 10.) == 10.
    @test isa(initial_partition_function()[1], Array{Float64,2})
end
