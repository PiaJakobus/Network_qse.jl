module Io
using CSV
using DelimitedFiles
using DataFrames
export read_table

function read_table()
    #f = CSV.read("tables/part_frdm.asc")
    f = readdlm("tables/part_frdm.asc", ' ', AbstractString, skipstart = 3, skipblanks=true, )
    #strings = f[1:5:end]
    A = f[2:5:end]
    return A
end

print(size(read_table()))
print(read_table())
print(6)
end
