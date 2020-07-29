module Network_qse

#"""
# https://docs.julialang.org/en/v1/manual/unicode-input/
#"""
__precompile__(false)


using Optim
using ForwardDiff
using Interpolations
using Plots
using NLsolve
using Dierckx

include("Io.jl")
include("Constants.jl")

export logsumexp
export find_nucl


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
how to write a test for this function?????
"""
function initial_partition_function(ω,A,Z,s,m)::Array{BigFloat,2}
    n_B = 1.0/const_m_B
    root_T⁻¹ = .√(1.0./(const_kmev.*data_T.*const_meverg))
    fp₀ = ω.*(2 .*s .+ 1)
    λ₀ = .√(const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_meverg/const_c^2)))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->(A[i]*fp₀[i]/(λ[i].^3.0))/n_B, 1:length(fp₀))...)
    return prefac
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
interpolation
"""
ω,A,Z,s,m =  Io.extract_partition_function()
pr = initial_partition_function(ω,A,Z,s,m)
npart = length(A)
inter_pr1 = [Spline1D(data_T, pr[j,:]) for j in 1:length(A)]
plot(LinRange(1e4,1e10,20),inter_pr1[end](LinRange(1e4,1e10,20)))
plot!(data_T, pr[end,:])


function log_charge_neutrality(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_kmev*T)
    prefact = abs.([inter_pr1[el](T) for el in 1:length(A)])
    result = log.(prefact.*(Z./A)./ ρ).+(μ[2] .* Z .+ μ[1] .* N .- E_b).*β
    return result
end

function log_mass_fraction(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    #μₙ,μₚ = μ
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_kmev*T)
    prefact = abs.([inter_pr1[el](T) for el in 1:length(A)])
    #println(prefact)
    result = log.(prefact./ ρ).+(μ[2] .* Z .+ μ[1] .* N .- E_b).*β
    return result
end



function mass_fraction(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    #μₙ,μₚ = μ
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_kmev*T)
    prefact = [inter_pr1[el](T) for el in 1:length(A)]
    result = ((prefact / ρ) .* exp.((μ[2] .* Z .+ μ[1] .* N .- E_b).*β))[:]
    return result
end


function logsumexp(arr)
    max = maximum(arr)
    dx = arr .- max
    sumexp = sum(exp.(dx))
    return max + log.(sumexp)
end



function f(F,x,T,yₑ,ρ,A,Z,m)
    F[1] = logsumexp(log_mass_fraction(x,T,ρ,A,Z,m))
    F[2] = logsumexp(log_charge_neutrality(x, T,ρ,A,Z,m))/log(yₑ) - 1
    return F
end



function ana_dev(J,μ, T,rho,y,A,Z,m)
    N = A .- Z
    J[1,1] = sum(N.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
    J[1,2] = sum(Z.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
    J[2,1] = sum((N.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y))
    J[2,2] = sum((Z.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y))
    return J
end

N = 10
sol_T = Array{BigFloat,2}(undef,(N,npart))
sol_T_auto = Array{BigFloat,2}(undef,(N,npart))
chempot = Array{BigFloat,2}(undef, (2,N))
chempot_auto = Array{BigFloat,2}(undef, (2,N))
eos_grid = eos((49.0,19.0,19.0))

rho = 1e8
y = 0.49
k = 0
t = 3e9

F = zeros(2)
test_dev1(x::Vector) = logsumexp(log_mass_fraction(x,t,rho,A,Z,m))
test_dev2(x::Vector) = logsumexp(log_charge_neutrality(x,t,rho,A,Z,m))/y-1.0
test_dev(x::Vector)  = f(F,x,t,y,rho,A,Z,m)
g = x -> [ForwardDiff.gradient(test_dev1,x),ForwardDiff.gradient(test_dev2,x)]
h = x -> ForwardDiff.jacobian(test_dev,x)
g(zeros(2))
h(zeros(2,2))

p = plot(g(zeros(2))[1])
for i in LinRange(1,10,100)
    plot!(p, g([i,i])[1])
end
p

for (i,t) in enumerate(LinRange(1e9,1e10,N)), (j,y) in enumerate(eos(range)[3][end])
    #sol = nlsolve((F,x)->f(F,x,t,y,rho), (J,x)->ana_dev(J,x,t,rho,y), [-0.2; 0.12])

    sol_auto = nlsolve((F,x)->f(F,x,t,y,rho,A,Z,m), [-0.2; 0.12],autodiff = :forward)#, method = :newton)#iterations = 1000)

    println("temp: ", t, sol_auto.zero, "  ", sum(mass_fraction([sol_auto.zero[1],sol_auto.zero[2]],t,rho,A,Z,m)),"  ", sum(exp.(log_charge_neutrality([sol_auto.zero[1],sol_auto.zero[2]],t,rho,A,Z,m)))-y)
    #sol_T[i,:] = mass_fraction([sol.zero[1],sol.zero[2]], t,rho)
    sol_T_auto[i,:] = mass_fraction([sol_auto.zero[1],sol_auto.zero[2]], t,rho,A,Z,m)
    #chempot[:,i] = sol.zero
    chempot_auto[:,i] = sol_auto.zero
end

plot(LinRange(1e6,1e10,N),sol_T .+ 0.0000001, yaxis=:log, ylims=(10e-6,1.0),xlabel = "T [K]", ylabel = "Xᵢ")
plot(LinRange(1e6,1e10,N),sol_T_auto .+ 0.0000001, yaxis=:log, ylims=(10e-6,1.0),xlabel = "T [K]", ylabel = "Xᵢ")

plot(LinRange(1e6,1e10,N),chempot[1,:], seriestype = :scatter, xlabel = "T [K]", ylabel ="μₙ", legend = :false)
plot!(LinRange(1e6,1e10,N),chempot_auto[1,:], seriestype = :scatter, xlabel = "T [K]", ylabel ="μₙ", legend = :false)

#seriestype = :scatter
plot(A,sol_T,seriestype = :scatter)
plot(G[3], transpose(sol_T),  xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)



function find_nucl(aa,zz,A,Z)
    y = findall(x->x == zz, Z)
    return y[findall(x->x==aa, A[y])]
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
function my_newton_raphson(μ,T,rho,A,Z,m)
    J = zeros(Float64, 2,2)
    y = 0.49
    F = Array{Float64,2}(undef, 2, 1)
    fun(x) = f(F,x,T,y,rho,A,Z,m)
    A, Z, m = G[[2,3,5]]
    N = A .- Z
    β = 1.0/(const_kmev*T)
    global ϵ = 1.0
    println("before loop ", typeof(μ))
    #global μⁱ⁺¹ = μ
    #μₙ,μₚ = μ
    global zaehler = 0
    while ϵ > 0.00001 && zaehler < 1000.0
        zaehler += 1
        println(">>> mu ", μ)
        #E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
        #prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
        J[1,1] = sum(N.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
        J[1,2] = sum(Z.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
        J[2,1] = sum((N.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y)) -1
        J[2,2] = sum((Z.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y)) -1
        det = J[1,1]*J[2,2] - J[1,2]*J[2,1]

        J⁻¹ = 1.0/det * [J[2,2] -J[2,1]; -J[1,2] J[1,1]]
        #println(">>> fun", fun(μ))
        #println(">>> J ", J)
        println(">>> det ", J)
        μⁱ⁺¹ = μ .- (zaehler / 1000.0).*[J[1,1]*fun(μ)[1] + J[1,2]*fun(μ)[2], J[2,1]*fun(μ)[1] + J[2,2]*fun(μ)[2]]
        μ = μⁱ⁺¹
        #println(">>> mu'", μ)
        #println()
        #print("-----",μⁱ⁺¹,"\n")
        ϵ = sum(mass_fraction([μⁱ⁺¹[1],μⁱ⁺¹[2]],T,rho))
        println(zaehler," ",">>> ϵ >>>", ϵ)
    end
    return μ
    end
end
soli = my_newton_raphson([0.1,0.2],3e9,2e8,A,Z,m)



end
