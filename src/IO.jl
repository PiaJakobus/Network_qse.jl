"""
read in follwing files:
----------------
mass-frdm95.dat:
----------------
Ground state properties
based on the FRDM model
Format
------
Each record of the file contains:

   Z    : charge number
   A    : mass number
   El   : element symbol
   fl   : flag corresponding to 0 if no experimental data available
                                1 for a mass excess recommended by
                                  Audi&Wapstra (1995)
                                2 for a measured mass from
                                  Audi&Wapstra (1995)
   Mexp : experimental or recommended atomic mass excess in MeV of
          Audi&Wapstra (1995)
   Mth  : calculated FRDM atomic mass excess in MeV
   Emic : calculated FRDM microscopic energy in MeV
   beta2: calculated quadrupole deformation of the nuclear ground-state
   beta3: calculated octupole deformation of the nuclear ground-state
   beta4: calculated hexadecapole deformation of the nuclear ground-state
   beta6: calculated hexacontatetrapole deformation of the nuclear
          ground-state

The corresponding FORTRAN format is (2i4,1x,a2,1x,i1,3f10.3,4f8.3)

----------------
part_frdm.asc:
----------------
Each of the two files starts with 4 header lines briefly
summarizing the contents of the given columns. This is followed by
entries sorted by charge and mass number of the isotope. Each
table ends with the line "END OF TABLE".
Each entry consists of 5 lines:
1. Isotope (in standard notation);
2. Charge number of isotope, mass number of isotope, ground state
   spin of the isotope;
3-5. Partition functions normalized to the g.s. spin;
   Third line: Partition functions for the  temperatures (in 10^9 K):
   0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7;
       Fourth line:  Partition functions for the temperatures (in 10^9 K):
   0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5;
       Fifth line: Partition functions for the temperatures (in 10^9 K):
   4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 .
Information for the next isotope starts after the last partition
function line.

>>>>>    PIA: ADDED partition functions and A,Z,s for n,p,d,T,3He
"""

function read_part_frdm()
    table_string = open("$(@__DIR__)/../tables/part_frdm.asc", "r") do f
        readlines(f)
    end
    str_f = split.(table_string, "\n")
    deleteat!(str_f, [1:4;])
    data_string = str_f[2:5:length(str_f)]

    k = vcat(map(x->split.(x, " "), data_string)...)
    k1 = map.(j-> filter!(i -> i != "", k[j]), 1:length(k))
    k2 = map(n->tryparse.(Float64,n), k1)
    k3 = permutedims(hcat(k2...))

    G_string = vcat(map(n-> str_f[3+5*n:5+5*n], 0:floor(Int, length(str_f)/5)-1)...)
    G_substring = vcat(map(x->split.(x, " "), G_string)...)
    G1 = map.(j -> filter!(i -> i != "", G_substring[j]), 1:length(G_substring))
    G2 = map(n->tryparse.(Float64,n), G1)
    G3 = permutedims(hcat(G2...))
    G4 = map(n->G3[1+3*n:3+3*n, 1:8], 0:floor(Int, size(G3)[1]/3)-1)
    G5 = map(x -> hcat(transpose(x)...), G4)
    return k3, G5
end


"""
    read_mass_frdm()
read ground state energies
"""
function read_mass_frdm()
    table_string = open("$(@__DIR__)/../tables/mass-frdm95.dat", "r") do f
        readlines(f)
    end
    b = split.(table_string, "\n")
    deleteat!(b, [1:4;])
    k = vcat(map(x->split.(x, " "), b)...)
    k1 = map.(j-> filter!(i -> i != "", k[j]), 1:length(k))
    k2 = map(i->deleteat!(i, 3), k1)
    k3 = map(n->tryparse.(Float64,n), k2)
    #k4 = map(x -> x[3] == 0.0 ? x[1:5] : x[1:4] , k3) # automatically parsed like this...
    k4 = map(x->x[1:4] , k3)
    k4_arr = permutedims(hcat(k4...))
    #k5 = k4_arr[k4_arr[:,3] .!= 0.0, :] # discard flag = 0
    return k4_arr
end



"""
    read_species()

reads species from table
"""
function read_species()
    strings = open("$(@__DIR__)/../tables/species.txt", "r") do f
        readlines(f)
    end
    number_species = parse(Int,strings[1])
    splitting = split.(strings, "\n")
    deleteat!(splitting, [1])
    k = vcat(map(x->split.(x, " "), splitting)...)
    k1 = map.(j-> filter!(i -> i != "", k[j]), 1:length(k))
    k2 = map(i->deleteat!(i, 1), k1)
    k3 = map(n->tryparse.(Float64,n), k2)
    k4 = permutedims(hcat(k3...))
    return k4, number_species
end


"""
    extract_partition_function()

checks indentical charge and atomic numbers
in mass_frdm and part_frdm.
Interpolates ω(T)
"""
function extract_partition_function()
    d1 = read_mass_frdm()
    d2, g = read_part_frdm()
    m_charge_number = Int64.(d1[:,1])
    m_atomic_number = Int64.(d1[:,2])
    m_mass          = d1[:,4]
    p_charge_number = Int64.(d2[:,1])
    p_atomic_number = Int64.(d2[:,2])
    T = 1e9.*Float64[0.01, 0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]
    #TODO: write a config where I store hard coded things
    result        = Array{AtomicProperties, 1}(undef, 0)
    for i in eachindex(m_charge_number)
        for j in eachindex(p_charge_number)
            if (((m_charge_number[i] == p_charge_number[j]) && (m_atomic_number[i] == p_atomic_number[j])))
                atomProp = AtomicProperties(m_charge_number[i],
                            m_atomic_number[i],
                            d2[j,3],
                            m_mass[i],
                            t-> Network_qse.LinearInterpolation(T, vcat(g[j]...))(t))
                push!(result, atomProp)
            #else # use all elements of part table, incl with no partition function, set it to 1!
            #    atomProp = AtomicProperties(m_charge_number[i],
            #                m_atomic_number[i],
            #                0.0,
            #                m_mass[i],
            #                t-> 1.0)
            #    push!(result, atomProp)
            end
        end

    end

    return result
end
