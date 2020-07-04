module Network_qse

__precompile__(false)

using CSV
using DelimitedFiles
using DataFrames

export read_part_frdm
export read_species
export read_mass_frdm
export extract_partition_function
export initial_partition_function
export findnearest
export linear_interpolation

include("Io.jl")


function testing()
    return 1
end

initial_partition_function()




end
