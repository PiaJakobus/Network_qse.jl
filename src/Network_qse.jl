module Network_qse
using CSV
using DelimitedFiles
using DataFrames
export mass
export initial_partition_function
include("Io.jl")

Nₐ = 1.
k = 1.
h = 1.
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
    n_B = ρ/const_m_B
    β   = 1.0/(const_kmev * T)
    scr3 = 1. # seems odd in code, not set..later sum(x_i)
    ω = G[1] # careful dimensions 3 x 8 matrix
    A = G[2]
    Z = G[3]
    s = G[4]
    m = G[5]
    N = A - Z
    scr3 = 1.0
    scr2 = findnearest(data_T, T)

    ω₁ = permutedims(hcat(map.(x -> ω[x][1, 1], [1:npart])...))
    ω₂ = permutedims(hcat(map.(x -> ω[x][1, 2], [1:npart])...))
    fp₀ = (ω₁ * scr2  .+ ω₂ * scr3).*(2.0*s .+ 1.0)
    scr1 = .√(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2)/const_hh.^2).^3
    binding_energy =  m - Z*m[2] + N*m[1]
    f₀ = fp₀ .* scr1 .* A / n_B
    return f₀
end

initial_partition_function(1,1)




end
