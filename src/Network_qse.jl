module Network_qse

"""
https://docs.julialang.org/en/v1/manual/unicode-input/
"""
__precompile__(false)

using CSV
using DelimitedFiles
using DataFrames
using Roots
using Optim
using ForwardDiff
using Interpolations
using Plots
using NLsolve
using Dierckx

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
export charge_neutrality
export mass_i
export mass_fraction

struct my_data
          A::Array{BigFloat}
          Z::Array{BigFloat}
          G::Vector
          eos::Vector{BigFloat}
          masses::Array{BigFloat}
          range::Tuple{BigFloat}
end

function eos(ind)::Vector
    i = 1.0:ind[1]
    j = 1.0:ind[2]
    k = 1.0:ind[3]
    rho = 10.0.^(log10(1e6) .+ i./49.0 .* log10(1e10/1e6))
    tem = 10.0.^(log10(2e9) .+ j./19.0 .* log10(9.9e9/2e9))
    y_e = collect(0.5 .+ (k.-19.0)./19.0 .* (0.5-0.405))
    return [rho,tem,y_e]
end

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
    ω,A,Z,s,m = G[1:5]
    root_T⁻¹ = .√(1.0./(const_kmev.*data_T.*const_meverg))
    fp₀ = ω.*(2 .*s .+ 1)
    λ₀ = .√(const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_meverg/const_c^2)))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->(A[i]*fp₀[i]/(λ[i].^3.0))/n_B, 1:length(fp₀))...)
    return prefac
end
initial_partition_function()





G_all = extract_partition_function()
range = (49.0,19.0,19.0)
eos_grid = eos((49.0,19.0,19.0))
m_n = 8.071
m_p = 7.289
#iron_iso = findall(x->x==56, G_all[2])
pr = initial_partition_function()#[iron_iso,:]
npart = length(G[1])
ω,A,Z,s,m = G[1:5]
plot(m,[G[1][j][1] for j in 1:length(Z)],seriestype = :scatter)


"""
interpolation
"""
inter_pr3 = [LinearInterpolation(data_T, pr[j,:]) for j in 1:length(G[2])]
inter_pr1 = [Spline1D(data_T, pr[j,:]) for j in 1:length(G[2])]
plot(data_T,inter_pr3[1000](data_T))
plot(data_T,inter_pr1[5000](data_T))
plot!(data_T, pr[10,:])


function charge_neutrality(μ,T::Float64,ρ::Float64)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    #μₙ,μₚ = μ
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_kmev*T)
    prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
    if any(x->x<=0, prefact)
        "domain error"
    end
    #println(prefact)
    #result = ((prefact.*(Z./A) / ρ) .* exp.((μₚ .* Z .+ μₙ .* N .- E_b).*β))[:]
    result = logsumexp(log.(prefact.*(Z./A)./ ρ).+(μ[2] .* Z .+ μ[1] .* N .- E_b).*β)
    return result
end

function mass_fraction(μ,T::Float64,ρ::Float64)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    #μₙ,μₚ = μ
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_kmev*T)
    prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
    #result = log.(prefact./ ρ).+ logsumexp.((μₚ .* Z .+ μₙ .* N .- E_b).*β)
    #println(prefact)
    if any(x->x<=0, prefact)
        "domain error"
    end
    result = logsumexp(log.(prefact./ ρ).+(μ[2] .* Z .+ μ[1] .* N .- E_b).*β)
    #result = (prefact./ ρ).*logsumexp((μₚ .* Z .+ μₙ .* N .- E_b).*β)
    return result
end

prefact = [inter_pr1[el](3e9) for el in 1:length(G[2])]


function mass_i(μ,T::Float64,ρ::Float64)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    #μₙ,μₚ = μ
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_kmev*T)
    prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
    result = ((prefact / ρ) .* exp.((μ[2] .* Z .+ μ[1] .* N .- E_b).*β))[:]
    return result
end



function logsumexp(arr)
    max = maximum(arr)
    dx = arr .- max
    sumexp = sum(exp.(dx))
    return max + log.(sumexp)
end

mass_i([1.1,1.1],1.0,2.0)


function f!(F,x,T,yₑ,ρ)
    F[1] = mass_fraction(x,T,ρ)
    F[2] = charge_neutrality(x, T,ρ)/log(yₑ) - 1
    return F
end



function ana_dev(J,μ, T,rho,y)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    J[1,1] = sum(N.*mass_i(μ, T,rho))/sum(mass_i(μ, T,rho))
    J[1,2] = sum(Z.*mass_i(μ, T,rho))/sum(mass_i(μ, T,rho))
    J[2,1] = sum((N.*Z./A).*mass_i(μ, T,rho))/(sum((Z./A).*mass_i(μ, T,rho))*log(y))
    J[2,2] = sum((Z.*Z./A).*mass_i(μ, T,rho))/(sum((Z./A).*mass_i(μ, T,rho))*log(y))
    return J
end


sol_T = Array{BigFloat,2}(undef,(20,npart))
chempot = Array{BigFloat,2}(undef, (2,20))
rho = 1e8
y = 0.49
k = 0
t = 3e9


sum(mass_i([1.5e-5,1.1e-5], 3e9, 1e9))
exp(mass_fraction([1.5e-5,1.1e-5], 3e9, 1e9))

data_T
data_T = 1e9.*Float64[0.01, 0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]


for (i,t) in enumerate(LinRange(1e9,1e10,10)), (j,y) in enumerate(eos(range)[3][end])
#for (i,t) in enumerate(LinRange(eos(range)[2][1],eos(range)[2][end],20)), (j,y) in enumerate(eos(range)[3][1])

    fun(F,x) = f!(F,x,t,y,rho)

    ana!(J,x) = ana_dev(J,x,t,rho,y)
    sol = nlsolve(fun, ana!, [0.01; 0.02])

    #jacobus!(J,x) = jacobian(J,x,t,rho,y)
    #sol = nlsolve(fun, jacobus!, [0.00001; 0.000012])

    #sol = nlsolve(fun, [0.00001; 0.000012],autodiff = :forward, method = :newton)#iterations = 1000)

    println("temp: ", t, sol.zero, "  ", exp(mass_fraction([sol.zero[1],sol.zero[2]],t,rho)) -1.0,"  ", exp(charge_neutrality([sol.zero[1],sol.zero[2]],t,rho))-0.49)
    sol_T[i,:] = mass_i([sol.zero[1],sol.zero[2]], t,rho)
    chempot[:,i] = sol.zero
end

#method = :anderson, method = :newton
#autodiff = :forward,
tresh = map(y->findall(x->x > 0.001, sol_T[y,:]), 1:20)

sol_T[:,iron_iso]
fe   = findall(x->x == 26, G_all[3])
chr  = findall(x->x == 24, G_all[3])
cob  = findall(x->x == 27, G_all[3])
ni   = findall(x->x == 28, G_all[3])
cop  = findall(x->x == 29, G_all[3])
ti   = findall(x->x == 22, G_all[3])

fe56  = fe[findall(x->x==56, A[fe])]
fe54  = fe[findall(x->x==50, A[findall(x->x == 26, G_all[3])])]
chr52 = chr[findall(x->x==52, A[findall(x->x == 24, G_all[3])])]
cob55 = cob[findall(x->x==55, A[findall(x->x == 27, G_all[3])])]
ni56  = ni[findall(x->x==56, A[findall(x->x == 28, G_all[3])])]
cop55 = cop[findall(x->x==55, A[findall(x->x == 29, G_all[3])])]
ti50  = ti[findall(x->x==50, A[findall(x->x == 22, G_all[3])])]

sol_T[1:10,:]

plot(LinRange(eos(range)[2][1],eos(range)[2][end],10),sol_T[1:10,:] .+ 0.0000001, yaxis=:log, ylims=(10e-3,1.0),xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)

plot(LinRange(eos(range)[2][1],eos(range)[2][end],10),sol_T[1:10,fe56] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)
plot!(LinRange(eos(range)[2][1],eos(range)[2][end],20),sol_T[:,fe54] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)
plot!(LinRange(eos(range)[2][1],eos(range)[2][end],20),sol_T[:,chr52] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)
plot!(LinRange(eos(range)[2][1],eos(range)[2][end],20),sol_T[:,cob55] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)
plot!(LinRange(eos(range)[2][1],eos(range)[2][end],20),sol_T[:,ni56] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)
plot!(LinRange(eos(range)[2][1],eos(range)[2][end],20),sol_T[:,cop55] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)
plot!(LinRange(eos(range)[2][1],eos(range)[2][end],20),sol_T[:,ti50] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)



plot(LinRange(eos(range)[2][1],eos(range)[2][end],200),chempot[1,:] + chempot[2,:], seriestype = :scatter, xlabel = "T [K]", ylabel ="μₙ", legend = :false)
#seriestype = :scatter
plot(G[3], transpose(sol_T),  xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)

sol_T[2]
chempot[1,:] + chempot[2,:]

eos(range)


"""
https://github.com/JuliaMath/Roots.jl
"""
function find_root()
    f(x::Vector) = teste(x,1.5e9,eos_grid[1][end])
    x = optimize(f, [1.5e-5,1.1e-5],Newton())
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
function my_newton_raphson(μ,T,rho)
    J = zeros(Float64, 2,2)
    y = 0.49
    F = Array{Float64,2}(undef, 2, 1)
    fun(x) = f!(F,x,T,y,rho)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    β = 1.0/(const_kmev*T)
    global ϵ = 1.0
    println("before loop ", typeof(μ))
    #global μⁱ⁺¹ = μ
    #μₙ,μₚ = μ
    while ϵ > 0.00001
        println(">>> mu ", μ)
        #E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
        #prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
        J[1,1] = sum(N.*mass_i(μ, T,rho))/sum(mass_i(μ, T,rho))
        J[1,2] = sum(Z.*mass_i(μ, T,rho))/sum(mass_i(μ, T,rho))
        J[2,1] = sum((N.*Z./A).*mass_i(μ, T,rho))/(sum((Z./A).*mass_i(μ, T,rho))*log(y))
        J[2,2] = sum((Z.*Z./A).*mass_i(μ, T,rho))/(sum((Z./A).*mass_i(μ, T,rho))*log(y))
        det = J[1,1]*J[2,2] - J[1,2]*J[2,1]

        J⁻¹ = 1.0/det * [J[2,2] -J[2,1]; -J[1,2] J[1,1]]
        println(">>> fun", fun(μ))
        println(">>> J ", J)
        println(">>> det ", det)
        #μⁱ⁺¹ = μ .- 0.01.*J⁻¹.*fun(μ)
        μⁱ⁺¹ = μ .- 0.01.*[J[1,1]*fun(μ)[1] + J[1,2]*fun(μ)[2],J[2,1]*fun(μ)[1] + J[2,2]*fun(μ)[2]]
        μ = μⁱ⁺¹
        println(">>> mu'", μ)
        println()
        #print("-----",μⁱ⁺¹,"\n")
        ϵ = exp(mass_fraction([μⁱ⁺¹[1],μⁱ⁺¹[2]],T,rho)) -1.0
        #println(μⁱ⁺¹)
    end
    return μ
    end
end
soli = my_newton_raphson([0.001,0.002],3e9,2e8)



end
