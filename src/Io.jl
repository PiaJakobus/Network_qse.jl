export read_part_frdm
export read_mass_frdm
export extract_partition_function
f = "tables/part_frdm.asc"
# https://docs.julialang.org/en/v1/manual/networking-and-streams/

export initial_partition_function
export findnearest
export linear_interpolation

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
data_T = Float64[0.01, 0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]

function initial_partition_function()
    """
            out: prefactor of X_i, using Boltzman stat.
            http://cococubed.asu.edu/code_pages/nse.shtml
            no interpolation, no ρ
                 prefac  = A*ω*λ⁻³
                 E_b = m - Z * mₚ - N * mₙ
            degeneracy gᵢ in ωᵢ = 2J(Zᵢ, Aᵢ, Eᵢ) + 1
    """

    n_B = 1.0/const_m_B
    β   = 1.0/(const_kmev * 1.0)
    scr3 = 1. # seems odd in code, not set..later sum(x_i)
    ω,A,Z,s,m = G[1:5] # careful dimensions G[1] isa 3 x 8 matrix
    N = A - Z
    root_T⁻¹ = .√(1.0./data_T)
    #i = findnearest(data_T, T)[1]
    #i_1 = Int8[ceil(i/8), i%8] # map i∈[1..n*m] to [n,m]
    fp₀ = map(x->ω[x]*(2.0*s[x]+1.0), [1:length(s);])
    λ₀ = .√const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2))
    λ = root_T⁻¹*λ₀
    E_b =  m - Z*m[2] + N*m[1]
    prefac = map(i->(A[i]*fp₀[i]*λ[i].^3.0)/n_B, [1:length(fp₀);])
    return prefac
end

test = initial_partition_function()

function parition_function()
    """
            non-interacting ideal gases,
            see Saha equation, using BS
            E_binding required as input, see
            Finite Range Droplet Model (FRDM)

            G(Z,A,T) = ∑ᵢ gᵢ(Zᵢ,Aᵢ) exp(-Eᵢ/kT)
    """
end

function findnearest(a,x)
       n = length(a)
       n > 0 || return 0:-1
       i1 = searchsortedlast(a,x)
       i0 = i1
       if i1>0
           while i0>1 && a[i0-1]==a[i0]
               i0 -= 1
           end
           d = x-a[i1]
       else
           i0 = i1+1
           d = a[i1+1]-x
       end
       i2 = i1
       if i2<n && a[i2+1]-x<d
           i0 = i2+1
           d = a[i2+1]-x
           i2 += 1
       end
       while i2<n && a[i2+1]-x==d
           i2 += 1
       end
       return i0:i2
end

function linear_interpolation(xₐᵣᵣ, yₐᵣᵣ, x)
    """
    forward interpolation
    xₐᵣᵣ = [xᵢ,xᵢ₊₁]
    yₐᵣᵣ = [yᵢ,yᵢ₊₁]
    """
    y⁺ = yₐᵣᵣ[1] + (x - xₐᵣᵣ[1])*(yₐᵣᵣ[2] - yₐᵣᵣ[1])/(xₐᵣᵣ[2] - xₐᵣᵣ[1])
    return y⁺
end




findnearest(data_T, 6)
