module Network_qse

__precompile__(false)

using CSV
using DelimitedFiles
using DataFrames
using Roots

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
    ω,A,Z,s,m = G[1:5]
    root_T⁻¹ = .√(1.0./data_T)
    #fp₀ = map(om->(2.0.*s.+1.0)*om, ω)
    fp₀ = ω.*(2 .*s .+ 1)
    λ₀ = .√const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->(A[i]*fp₀[i]*λ[i].^3.0)/n_B, 1:length(fp₀))...)
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



function saha_equation(μₙ,μₚ,T,ρ)
    prefactor = initial_partition_function()
    A, Z, m = G[[2,3,5]]
    N = A - Z
    E_b =  m - Z*m[2] + N*m[1]
    xᵢ = exp.(μₚ .* z)
    return [x^2, x^3]
end

μₚ = 1
μₙ = 1
T = 1
A, Z, m = G[[2,3,5]]
N = A - Z
E_b =  m - Z*m[2] + N*m[1]
pr = initial_partition_function()
xᵢ = pr .* exp.((μₚ .* Z + μₙ .* N - E_b).*(const_k_B * T))
Yₑ = xᵢ .* Z ./ A

function find_root()
    f(μ) = sum(saha_equation(μₙ,μₚ,T,ρ)) - 1
    g(μ) = sum(charge_neutrality())
    find_zeros(f,-∞, ∞)
end



initial_partition_function()




end
