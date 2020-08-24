using Test
using Network_qse


@testset "network_nse" begin
    include("Network_qse.jl")
end

@testset "tools" begin
    include("Tools.jl")
end

@testset "io" begin
    include("IO.jl")
end

@testset "DataTypes" begin
    include("DataTypes.jl")
end

@testset "Boltzmann" begin
    include("Boltzmann.jl")
end
