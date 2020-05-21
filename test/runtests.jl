using Test
#using Network_qse
#include("../src/Io.jl")

@testset "first_try" begin
    @test 1==1
    @test foo(1) == 1
end

@testset "io" begin
    @test isa(read_table(), Array{String})
end
