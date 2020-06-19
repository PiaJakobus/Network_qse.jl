using Test
using Network_qse
using DelimitedFiles

include("../src/Io.jl")

@testset "first_try" begin
    @test 1==1
    @test foo(1) == 1
end

@testset "io" begin
    @test isa(read_part_frdm(), Tuple{Array{Array{Float64,1},1},Array{Array{Float64,1},1}})
    @test isa(read_mass_frdm(), Array{Array{Float64,1},1})
end
