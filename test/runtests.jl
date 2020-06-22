using Test
#using Network_qse
using DelimitedFiles

include("../src/Io.jl")

@testset "first_try" begin
    @test 1==1

    @test mass() == 1
end

@testset "io" begin
    @test isa(read_part_frdm(), Tuple{Array{Float64,2},Array{Array{Array{Float64,1},1},1}})
    @test isa(read_mass_frdm(), Array{Float64,2})
    @test isa(extract_partition_function(), Tuple{Array{Float64,2},Array{Float64,1},Array{Float64,1}})
end
