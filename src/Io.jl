export read_part_frdm
export read_mass_frdm
export extract_partition_function
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
    part_arr = reshape(hcat(G_mod...), (length(G_mod[1]), length(G_mod)))
    return data_arr, part_arr
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

    return k4_arr
end

d1 = read_mass_frdm()
d2, g = read_part_frdm()
m_charge_number = d1[:,1]
m_atomic_number = d1[:,2]
p_charge_number = d2[:,1]
p_atomic_number = d2[:,2]
m_zz_aa = d1[:,[1,2]]
p_zz_aa = d2[:,[1,2]]


function extract_partition_function()
    fpart = Array{Array{Float64,1},1}[]
    atomic_number = Vector{Float64}()
    charge_number = Vector{Float64}()
    for i in eachindex(m_charge_number)
        for j in eachindex(p_charge_number)
            if (m_charge_number[i] == p_charge_number[j]) && (m_atomic_number[i] == p_atomic_number[j])
                push!(fpart, g[:,j])
                push!(atomic_number, m_atomic_number[i])
                push!(charge_number, m_charge_number[i])
            end
        end
    end
    #f_part = reshape(hcat(fpart[1]...), (length(fpart[1][1]), length(fpart[1])))
    return fpart, atomic_number, charge_number
end




#A_in_p = in.(i,p_charge_number)
