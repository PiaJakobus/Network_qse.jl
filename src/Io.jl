export read_part_frdm
export read_mass_frdm
export extract_partition_function
f = "tables/part_frdm.asc"
# https://docs.julialang.org/en/v1/manual/networking-and-streams/

export initial_partition_function


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
m_mass          = d1[:,3]
m_spin          = d1[:,4]
p_charge_number = d2[:,1]
p_atomic_number = d2[:,2]
m_zz_aa = d1[:,[1,2]]
p_zz_aa = d2[:,[1,2]]

function read_species()
    string = open("tables/species.txt", "r") do f
        readlines(f)
    end
    number_species = parse(Int,string[1])
    splitting = split.(string, "\n")
    deleteat!(splitting, [1])
    k = vcat(map(x->split.(splitting[x], " "), [1:length(splitting);])...)
    k1 = map(n->tryparse.(Float64,k[n]), [1:length(k);])
    k2 = map.(i-> filter!(x->x≠nothing,k1[i]), [1:length(k1);])
    k3 = map(x->identity.(k2[x]), [1:length(k2);])
    k4 = permutedims(hcat(k3...))
    return k4, number_species
end

function extract_partition_function()
    fpart         = Array{Float64,2}[]
    atomic_number = Vector{Float64}()
    charge_number = Vector{Float64}()
    spin          = Vector{Float64}()
    mass          = Vector{Float64}()
    for i in eachindex(m_charge_number)
        for j in eachindex(p_charge_number)
            if (m_charge_number[i] == p_charge_number[j]) && (m_atomic_number[i] == p_atomic_number[j])
                push!(fpart, transpose(reshape(hcat(g[:,j]...), (length(g[:,j][1]), length(g[:,1])))))
                push!(atomic_number, m_atomic_number[i])
                push!(charge_number, m_charge_number[i])
                push!(spin, m_spin[i])
                push!(mass, m_mass[i])
            end
        end
    end
    return fpart, atomic_number, charge_number, spin, mass
end



G = extract_partition_function()
npart = length(G[1])
masses = read_mass_frdm()
nmass = length(masses[:,1])
const_m_B = 1.66e-24 # baryon mass
const_kmev = 8.61829e-11
const_meverg = 1.602e-6
const_k_B = 1.380658e-16
const_c = 2.99792458e10
const_h_barc = 197.327e-13
const_hh = const_h_barc / const_c * 2.0 * π * const_meverg

function initial_partition_function(ρ,T)
    n_B = 1.0/const_m_B
    β   = 1.0/(const_kmev * 1.0)
    scr3 = 1. # seems odd in code, not set..later sum(x_i)
    ω = G[1] # careful dimensions 3 x 8 matrix
    A = G[2]
    Z = G[3]
    s = G[4]
    m = G[5]
    N = A - Z
    scr3 = 1.0
    scr2 = 1.0

    ω₁ = permutedims(hcat(map.(x -> ω[x][1, 1], [1:npart])...))
    ω₂ = permutedims(hcat(map.(x -> ω[x][1, 2], [1:npart])...))
    fp₀ = (ω₁ * scr2  .+ ω₂ * scr3) .* transpose(2.0.*s .+ 1.0)
    scr1 = .√(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2)/const_hh.^2).^3
    binding_energy =  m - Z*m[2] + N*m[1]
    f₀ = fp₀ .* transpose(scr1) .* transpose(A) / n_B
    return vcat(f₀...)
end

test = initial_partition_function(1,1)
