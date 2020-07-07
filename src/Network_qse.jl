module Network_qse

__precompile__(false)

using CSV
using DelimitedFiles
using DataFrames
using Roots
using Optim
using ForwardDiff


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



"""
        initial_partition_function()
returns the prefactor of X_i, as
formulated within Boltzman statistics
in Saha equation, see i.g. here
http://cococubed.asu.edu/code_pages/nse.shtml
https://en.wikipedia.org/wiki/Saha_ionization_equation
prefac = A*ω*(2s+1)/λ³
-- M   = A*m_B + m / c²
-- mB  = mₚ = mₙ, m = ground state mass of nuclei
-- λ   ≡ √h²/2πMβ
-- fp₀ = ω*g
-- g   = 2s + 1
Don't forget the ρ !
"""
function initial_partition_function()
    n_B = 1.0/const_m_B
    β   = 1.0/(const_kmev * 1.0)
    ω,A,Z,s,m = G[1:5]
    root_T⁻¹ = .√(1.0./data_T)
    fp₀ = ω.*(2 .*s .+ 1)
    λ₀ = .√const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->(A[i]*fp₀[i]/λ[i].^3.0)/n_B, 1:length(fp₀))...)
    return prefac
end



"""
call this only in main ... nothing happens here
except matrix multiplication from file inpput
later used in interpol(Float: T)
"""
pr = initial_partition_function()


"""
    interpol(T)
given a parameter T, find
corresponding value in multi-array
[[p₁¹,p₂¹,..,T¹],...,[p₁ⁱ,p₂ⁱ,..,Tⁱ]..]
interpolating between Tᵢ, Tᵢ₊₁
return pⁱⁿᵗᵉʳ. Contains list search
and linear interpolation scheme
i in our case: temperature grid (=24)
"""
function interpol(T::Float64)
    #searchsortedlast
    #i = findnearest(data_T, T)[1]
    #i_1 = Int8[ceil(i/8), i%8] # map i∈[1..n*m] to [n,m]
    return T
end



function saha_equation(μ::Array{Float64},T::Float64,ρ::Float64)
    A, Z, m = G[[2,3,5]]
    N = A - Z
    μₙ,μₚ = μ
    E_b =  m - Z*m[2] + N*m[1]
    prefact = interpol(T)
    xᵢ = prefact/ρ .* exp.((μₚ .* Z + μₙ .* N - E_b).*(const_k_B * T))
    yₑ = xᵢ .* Z ./ A
    return [xᵢ,yₑ]
end


f(x::Vector) = sum(sin, x)
x = rand(5)
g = x -> ForwardDiff.hessian(f, x); # g = ∇f
g(x)




function find_root(y::Float64, T::Float64, ρ::Float64)
    f(μ::Array{Float64}) = ∑(saha_equation(μ, T, ρ)) - [1, y]
    μ = my_newton_raphson(f, -10.0, 10.0)
end







end
