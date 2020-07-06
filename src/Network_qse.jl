module Network_qse

__precompile__(false)

using CSV
using DelimitedFiles
using DataFrames

include("Io.jl")
include("Constants.jl")
include("Tools.jl")

export read_part_frdm
export read_species
export read_mass_frdm
export extract_partition_function
export initial_partition_function
export parition_function
export findnearest
export linear_interpolation



G = extract_partition_function()
npart = length(G[1])
masses = read_mass_frdm()
nmass = length(masses[:,1])




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
    fp₀ = map(om->(2.0.*s.+1.0)*om, ω)
    λ₀ = .√const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2))
    λ = root_T⁻¹*λ₀
    E_b =  m - Z*m[2] + N*m[1]
    prefac = map(i->(A[i]*fp₀[i]*λ[i].^3.0)/n_B, 1:length(fp₀))
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
    #searchsortedlast
    #i = findnearest(data_T, T)[1]
    #i_1 = Int8[ceil(i/8), i%8] # map i∈[1..n*m] to [n,m]
    return 1
end






initial_partition_function()




end
