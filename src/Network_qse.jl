module Network_qse

"""
https://docs.julialang.org/en/v1/manual/unicode-input/
"""
__precompile__(false)

using CSV
using DelimitedFiles=
using DataFrames
using Roots
using Optim
using ForwardDiff
using Interpolations
using Plots


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
export my_newton_raphson
export saha_equation
export eos

struct my_data
          A::Array{Float64}
          Z::Array{Float64}
          G::Vector
          eos::Vector{Float64}
          masses::Array{Float64}
          range::Tuple{Float64}
end

G = extract_partition_function()
range = (49.0,19.0,19.0)
eos_grid = eos((49.0,19.0,19.0))
pr = initial_partition_function()
npart = length(G[1])
masses = read_mass_frdm()
nmass = length(masses[:,1])

inter_pr1 = [LinearInterpolation(data_T, pr[j,:]) for j in 1:length(G[2])]
inter_pr1 = mapslices.(x->LinearInterpolation(data_T, x), pr; dims = [1])

nodes = (data_T,)
inter_pr2 = [interpolate(nodes, pr[j,:], Gridded(Linear())) for j in length(G[2])]
plot(data_T, pr[4500,:])
for i=4500:10:4600
    p=display(plot!(data_T,pr[i,:]))
end
display(p)

"""
        initial_partition_function()
returns prefactor of X_i, as
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
    λ₀ = .√(const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_kmev/const_c^2)))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->(A[i]*fp₀[i]/λ[i].^3.0)/n_B, 1:length(fp₀))...)
    return prefac
end
initial_partition_function()


prefact = [inter_pr1[el](1.5e9) for el in 1:length(G[2])]


eos_grid[1]


function charge_neutrality(μ::Vector,T::Float64,ρ::Float64)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    μₙ,μₚ = μ
    E_b =  m .- Z*m[2] .+ N*m[1]
    β = 1.0/(const_kmev*T)
    prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
    result = (((prefact).*(Z./A) / (ρ)) .* exp.((μₚ .* Z + μₙ .* N - E_b).*(β)))[:]
    #result = prefact
    return result
end

y1 = charge_neutrality([1.1,2.1],1.5e9,eos_grid[1][end])
y2 = saha_equation([1.1e12,2.12],1.5e9,eos_grid[1][end])

plot(Z,y2)
plot(eos_grid[3], eos_grid[2])


function saha_equation(μ::Vector,T::Float64,ρ::Float64)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    μₙ,μₚ = μ
    E_b =  m .- Z*m[2] .+ N*m[1]
    β = 1.0/(const_kmev*T)
    prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
    result = ((prefact ./ ρ) .* exp.((μₚ .* Z + μₙ .* N - E_b).*(β)))[:]
    #yₑ = xᵢ .* Z ./ A
    return result
end





function eos(ind::Tuple{Float64})::Vector
    i = 1.0:ind[1]
    j = 1.0:ind[2]
    k = 1.0:ind[3]
    rho = 10.0.^(log10(1e6) .+ i./49.0 .* log10(1e10/1e6))
    tem = 10.0.^(log10(2e9) .+ j./19.0 .* log10(9.9e9/2e9))
    y_e = collect(0.5 .+ (k.-19.0)./19.0 .* (0.5-0.405))
    return [rho,tem,y_e]
end
eos(range)




"""
https://github.com/JuliaMath/Roots.jl
"""
function find_root()::Float64
    f(x::Float64) = x^2 - 1
    x = find_zeros(f, -10.0, 10.0)
end

"""
    Multivariate Newton raphson()
[xⁱ⁺¹₁..xⁱ⁺¹ₙ] = [xⁱ₁..xⁱₙ] - J⁻¹[f¹(xⁱ₁)..fⁿ(xⁱₙ)]
specificly:
calculate Hessian 2x2
[df[1]/dmun df[1]/dmup; df[2]/dumun df[2]/dmup]
J^-1 = 1/(ad-bc) * [d -b; -c a]
[mun,mup]' = [mun,mup] - J^-1 * f(mun,mup)
dXdμₙ   dXdμₚ
dYₑdμₙ  dYₑdμₚ
"""
function my_newton_raphson(μ::Array{Float64},T::Float64,ρ::Float64)::Array{Float64}
    y = 0.49
    G = extract_partition_function()
    f(μ::Array{Float64}) = collect(sum.([saha_equation(μ, T, ρ),charge_neutrality(μ,T,ρ)])) - [1, y]
    #G = extract_partition_function()
    A, Z = G[2:3]
    β = 1.0/(const_kmev*T)
    ϵ = 1.0
    μⁱ⁺¹ = [1.0,2.0]
    while ϵ > 0.01
        dXdμₙ  = sum((A .- Z)*β.*saha_equation(μ, T, ρ))
        dXdμₚ  = sum(β*Z.*saha_equation(μ, T, ρ))
        dYₑdμₙ = sum(dXdμₙ./(A.*Z))
        dYₑdμₚ = sum(dXdμₚ./(A.*Z))
        det = (dXdμₙ*dYₑdμₚ - dXdμₚ*dYₑdμₙ)
        J⁻¹ = 1.0/det * [dXdμₙ dXdμₚ; dYₑdμₙ dYₑdμₚ]
        μⁱ⁺¹ = μ .- J⁻¹.*f(μ)
        ϵ = √(abs(sum((μⁱ⁺¹ - μⁱ⁺¹).*(μⁱ⁺¹ - μⁱ⁺¹))))
    end
    return μⁱ⁺¹
    end
end
sol = my_newton_raphson([1.1,2.1],2.1e9,2.2e6)
saha_equation(sol[1,:],2.1e9,2.2e6)

end
