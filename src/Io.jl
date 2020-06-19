export read_part_frdm
export read_mass_frdm
f = "tables/part_frdm.asc"
# https://docs.julialang.org/en/v1/manual/networking-and-streams/

function read_part_frdm()
    table_string = open("tables/part_frdm.asc", "r") do f
        readlines(f)
    end
    str_f = split.(table_string, "\n")
    deleteat!(str_f, [1:4;])
    data_string = str_f[2:5:length(str_f)]
    data_substring = map(x->split.(data_string[x], " "), [1:length(data_string);])
    data_union = map(i->map(n->tryparse.(Float64,data_substring[i][n]), [1:length(data_substring[1]);]), [1:length(data_substring);])
    data_union1 = map.(i-> filter!(k->k≠nothing,data_union[i][1]), [1:length(data_union);])
    data_res = map(x->identity.(data_union1[x]), [1:length(data_union1);])

    G_string = vcat(map(n-> str_f[3+5*n:5+5*n], [0:floor(Int, length(str_f)/5)-1;])...)
    G_substring = vcat(map(x->split.(vcat(G_string[x]), " "), [1:length(G_string);])...)
    G_union = map(n->tryparse.(Float64,G_substring[n]), [1:length(G_substring);])
    G_union1 = map(y->filter!(x->x≠nothing,G_union[y]), [1:length(G_union);])
    G_res = map(x->identity.(G_union1[x]), [1:length(G_union1);])
    G_mod = map(n->G_res[1+3*n:3+3*n], [0:floor(Int, length(G_res)/3)-1;])

    data_arr = permutedims(reshape(hcat(data_res...), (length(data_res[1]), length(data_res))))
    part_arr = permutedims(reshape(hcat(G_mod...), (length(G_mod[1]), length(G_mod))))
    return data_res, G_mod
end


function read_mass_frdm()
    table_string = open("tables/mass-frdm95.dat", "r") do f
        readlines(f)
    end
    b = split.(table_string, "\n")
    deleteat!(b, [1:4;])
    k = vcat(map(x->split.(b[x], " "), [1:length(b);])...)
    k1 = map(n->tryparse.(Float64,k[n]), [1:length(k);])
    k2 = map.(i-> filter!(x->x≠nothing,k1[i]), [1:length(k1);])
    k3 = map(x->identity.(k2[x]), [1:length(k2);])
    k4 = map(i->k3[i][1:4], [1:length(k3);])
    k4_arr = permutedims(reshape(hcat(k4...), (length(k4[1]), length(k4))))

    return k4
end
