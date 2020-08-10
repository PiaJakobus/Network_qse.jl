module Network_qse

#"""
# https://docs.julialang.org/en/v1/manual/unicode-input/
#"""
__precompile__(false)


using Optim
using ForwardDiff
using Interpolations
using Plots
gr()
using NLsolve
using LsqFit
using Dierckx
using LinearAlgebra

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
    root_T⁻¹ = .√(1.0./(data_T))
    fp₀ = ω.*(2 .*s .+ 1)
    λ₀ = .√(const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_meverg/const_c^2)))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->((A[i]/ 6.02214076e23)*fp₀[i]/(λ[i].^3.0))/n_B, 1:length(fp₀))...)
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
inter_pr1 = [Spline1D(data_T, log.(pr[j,:])) for j in 1:npart]
plot(LinRange(1e7,1e10,50),exp.(inter_pr1[end](LinRange(1e7,1e10,50))),legend = :false)
plot!(data_T, pr[end,:],seriestype=:scatter,legend = :false)



function log_charge_neutrality(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(const_k_B*T)
    prefact = abs.([exp.(inter_pr1[el](T)) for el in 1:length(A)])
    result = log.(prefact.*(Z./A)./ ρ).+(μ[2] .* Z .+ μ[1] .* N .- E_b).*β
    return result
end

function log_mass_fraction(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    E_b =  (m .- Z*m_p .+ N*m_n)*const_meverg
    β = 1.0/(const_k_B*T)
    prefact = [exp.(inter_pr1[el](T)) for el in 1:length(A)]
    result = log.(prefact./ ρ).+(μ[2] .* Z .+ μ[1] .* N .- E_b).*β #μₙ,μₚ = μ
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

Z

function ana_dev(J,μ, T,rho,y,A,Z,m)
    N = A .- Z
    β = 1.0/(const_k_B*T)
    J[1,1] = sum(β.*N.*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/sum(exp.(log_mass_fraction(μ, T,rho,A,Z,m)))
    J[1,2] = sum(β.*Z.*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/sum(exp.(log_mass_fraction(μ, T,rho,A,Z,m)))
    J[2,1] = sum((β.*N.*Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/(sum((Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))*log(y))
    J[2,2] = sum((β.*Z.*Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/(sum((Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))*log(y))
    return J
end


F = zeros(2)
J = zeros(Float64, 2,2)
test_dev1(x::Vector) = logsumexp(log_mass_fraction(x,t,rho,A,Z,m))
test_dev2(x::Vector) = logsumexp(log_charge_neutrality(x,t,rho,A,Z,m))/log(y) - 1
test_dev(x::Vector)  = f(F,x,t,y,rho,A,Z,m)
g = x -> [ForwardDiff.gradient(test_dev1,x),ForwardDiff.gradient(test_dev2,x)]
h = x -> ForwardDiff.jacobian(test_dev,x)
g(zeros(2))
h(zeros(2,2))

p = plot(g([1000,1000])[2])
for i in LinRange(0,1000,5)
    plot!(p, g([i,i])[1],color=:red, legend = :false)
    plot!(p, g([i,i])[2],color=:green, legend = :false)
    plot!(p, ana_dev(J,[i,i],t,rho,y,A,Z,m)[1,:],color=:orange,seriestype = :scatter,legend = :false)
    plot!(p, ana_dev(J,[i,i],t,rho,y,A,Z,m)[2,:],color=:blue,seriestype = :scatter,legend = :false)
end
p



eos_grid     = eos((49.0,19.0,19.0))
x_guess = [6.492724770496382e-6, -8.232974345747265e-6]

function solve(x_guess, f_root,dev_a)
    N            = 1
    N_y          = 40
    N_rho        = 1
    sol_T        = Array{BigFloat,4}(undef,(N,N_y,N_rho,npart))
    sol_T_NR       = Array{BigFloat,4}(undef,(N,N_y,N_rho,npart))
    sol_T_auto   = Array{BigFloat,4}(undef,(N,N_y,N_rho,npart))
    chempot      = Array{BigFloat,4}(undef, (2,N,N_y,N_rho))
    chempot_auto = Array{BigFloat,4}(undef, (2,N,N_y,N_rho))
    chempot_NR   = Array{BigFloat,4}(undef, (2,N,N_y,N_rho))
    t = 9e9
    rho = 1e7
    y = 0.41
    i = 1
    j = 1
    k = 1
    for (j,y) in enumerate(LinRange(0.41,0.6,N_y))#, (k,rho) in enumerate(LinRange(1e9,1e9,N_rho))
    #for (i,t) in enumerate(LinRange(1e9,12e10,N))
        sol_auto = nlsolve((F,x)->f_root(F,x,t,y,rho,A,Z,m),x_guess,autodiff = :forward)#, method = :newton)#iterations = 1000)
        #sol_NR   = my_newton_raphson(x,t,y,rho,A,Z,m)
        #sol = nlsolve((F,x)->f(F,x,t,y,rho,A,Z,m), (J,x)->ana_dev(J,x,t,rho,y,A,Z,m), [6.492724770496382e-6, -8.232974345747265e-6])
        println("chempot_auto: ", sol_auto.zero)#, "chempot_NR: ", sol_NR)
        #println("chempot_NR: ", sol_NR, f(F,sol_NR,t,y,rho,A,Z,m))
        #println("chempot_ana: ", sol.zero)
        println("temp: ", t, " << >> ","y_e ",y," ", sum(exp.(log_mass_fraction([sol_auto.zero[1],sol_auto.zero[2]],t,rho,A,Z,m)))-1,"  ", logsumexp(log_charge_neutrality([sol_auto.zero[1],sol_auto.zero[2]],t,rho,A,Z,m))/log(y) -1)
        #sol_T[i,j,k,:]      = exp.(log_mass_fraction([sol.zero[1],sol.zero[2]], t,rho,A,Z,m))
        #sol_T_NR[i,j,k,:]     = mass_fraction(sol_NR, t,rho,A,Z,m)
        sol_T_auto[i,j,k,:] = exp.(log_mass_fraction([sol_auto.zero[1],sol_auto.zero[2]], t,rho,A,Z,m))
        #chempot[:,i,j,k] = sol.zero
        chempot_auto[:,i,j,k] = sol_auto.zero
        x_guess = sol_auto.zero
        #chempot_NR[:,i,j,k] = sol_NR
    #end
    end
    return sol_T_auto,chempot_auto
end

X_i,mu = solve(x_guess,f,ana_dev)

X_i
function mapping()
    for (i,j) in enumerate(1:10;), (k,l) in enumerate(20:30;)
        println(j,l)
    end
end
mapping()

filter(x-> !=(x,map(j -> any(X_i[1,:,1,j] .> 0.001) && (X_i[1,:,1,,j), 1:length(Z))

filter(x -> !=(x,Bool), [1,2,3])


n > 0 || return 0:-1

range_T = LinRange(data_T[1],data_T[end],30)
fig_auto = plot(range_T,sol_T_auto[:,1,1,:] .+ 0.0000001,  ylims=(1e-2,1.0),yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]),legend = :false)
plot!(range_T,sol_T_NR[:,1,1,:] .+ 0.0000001, ylims=(1e-2,1.0), yaxis=:log,seriestype=:scatter,xlabel = "T [K]", ylabel = "Xᵢ",yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]),legend = :false, label=:"NR")
plot(range_T,sol_T[:,1,1,:] .+ 0.0000001, ylims=(1e-2,1.0), yaxis=:log,seriestype=:scatter,xlabel = "T [K]", ylabel = "Xᵢ",yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]),legend = :false, label=:"ana")

savefig(fig_auto,"nse_auto_NR_1e7.pdf")
savefig(fig_auto,"nse_auto_analyt.pdf")


dict2 = Dict("n" => [1,0], "p" => [1,1], "he3" => [3,2], "o16" => [16,8], "ne20" => [20,10], "mg24" => [24,12],
"si28" => [28,14], "s32" => [32,16], "ar36" => [36,18], "ca40" => [40,20], "ti44" => [44,22], "ti50" => [50,22],"cr48" => [48,24],
"chr52" => [52,24], "fe52" => [52,26], "fe54" => [54,26], "cob55" => [55,27], "cop55" => [55,29],"ni56" => [56,28], "fe56" => [56,26],
"v52" => [52,23], "ni62" => [62,28])

dict = Dict("n" => [1,0], "p" => [1,1], "he3" => [3,2], "o16" => [16,8], "ne20" => [20,10],
"si28" => [28,14], "ti50" => [50,22],"chr52" => [52,24], "fe54" => [54,26], "cob55" => [55,27], "cop55" => [55,29],"ni56" => [56,28], "fe56" => [56,26],
"ni62" => [62,28])




str = "ni62"
vary_T = plot(range_T,sol_T_auto[:,1,1,find_nucl(1,0,A,Z)[1]] .+ sol_T_auto[:,1,1,find_nucl(1,1,A,Z)[1]] .+ 0.0000001,yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = "p+n")#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
elem = plot(range_T,X_i[:,1,1,find_nucl(52,26,A,Z)[1]] .+ 0.0000001,xlabel = "T [K]", ylabel = "Xᵢ",label = str)
for k in keys(dict)
    plot!(range_T,sol_T_auto[:,1,1,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)
    #plot!(range_T,sol_T_NR[:,1,1,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001, seriestype =:scatter,yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
end
vary_T
savefig(elem,"ti50_den1e7.pdf")


fig_auto = plot(LinRange(0.41,0.5,N_y),sol_T_auto[1,:,1,:] .+ 0.0000001,  ylims=(1e-2,1.0),yaxis=:log,xlabel = "Yₑ", ylabel = "Xᵢ",yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]),legend = :false)
vary_y = plot(LinRange(0.41,0.6,40),X_i[1,:,1,find_nucl(52,26,A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "Yₑ", ylabel = "Xᵢ",label = "fe52")#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
for k in keys(dict)
    plot!(LinRange(0.41,0.6,40),X_i[1,:,1,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001,yaxis=:log, legend=:left,xlabel = "Yₑ", ylabel = "Xᵢ",label = k)
    #plot!(LinRange(0.41,0.5,N_y),sol_T[1,:,1,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001, seriestype =:scatter,yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
end
vary_y
savefig(fig_auto, "vary_ye_all_rho1e7_t3.5e8.pdf")
savefig(vary_y,"vary_ye_t3.5e9_rho1e9_higher_resol.pdf")




vary_rho = plot(LinRange(0.41,0.5,N_y),sol_T_auto[:,1,1,find_nucl(1,0,A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = "p")#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
for k in keys(dict)
    plot!(LinRange(0.41,0.5,N_rho),sol_T_auto[1,1,:,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)
    plot!(LinRange(0.41,0.5,N_rho),sol_T[1,1,:,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001, seriestype =:scatter,yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
end
vary_rho




plot(LinRange(1e8,1e9,N),sol_T[:,1,find_nucl(3,2,A,Z)[1]] .+ 0.0000001,color=:red, seriestype = :scatter, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ")
plot(LinRange(1e8,1e9,N),sol_T_auto[:,1,find_nucl(56,28,A,Z)[1]] .+ 0.0000001, color=:orange, yaxis=:log, xlabel = "T [K]", ylabel = "Xᵢ",title="3He")
savefig("he3.pdf")
#ylims=(10e-6,1.0),
plot(LinRange(1e8,1e9,N),chempot[1,:,1], seriestype = :scatter,color=:blue, xlabel = "T [K]", ylabel ="μₙ", legend = :false)
plot!(LinRange(1e8,1e9,N),chempot_auto[1,:,1], color=:green, xlabel = "T [K]", ylabel ="μₙ", legend = :false)
savefig("mu_n.pdf")

plot(vary_T,chempot_auto[2,:,1,1], seriestype = :scatter,color=:blue, xlabel = "T [K]", ylabel ="μₚ", legend = :false)
plot!(LinRange(1e8,1e9,N),chempot_auto[2,:,1], color=:green, xlabel = "T [K]", ylabel ="μₚ", legend = :false)
savefig("mu_p.pdf")


#seriestype = :scatter
plot(A,sol_T,seriestype = :scatter)
plot(G[3], transpose(sol_T),  xlabel = "T [K]", ylabel = "Xᵢ",legend = :false)


"""return index in (A/Z/m/s)"""
function find_nucl(aa,zz,A,Z)
    y = findall(x->x == zz, Z)
    return y[findall(x->x==aa, A[y])]
end

chempot_NR
sol_T_NR




range_T
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
soli = my_newton_raphson([7.862587588204352e-6, -9.617299012693205e-6],1e9,1e9,A,Z,m)
mat = ones(Float64, 2,2)
pinv(mat)
function my_newton_raphson(μ,T,y,rho,A,Z,m)
    mat = zeros(Float64, 2,2)
    #println(">>>><<<<", det(mat))
    #y = 0.41
    F = Array{Float64,2}(undef, 2, 1)
    fun(x) = f(F,x,T,y,rho,A,Z,m)
    N = A .- Z
    β = 1.0/(const_k_B*T)
    global ϵ = 1.0
    #global μⁱ⁺¹ = μ
    #μₙ,μₚ = μ
    global zaehler = 0
    while ϵ > 1e-10
        zaehler += 1
        #println(">>> mu ", μ)
        #E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
        #prefact = [inter_pr1[el](T) for el in 1:length(G[2])]
        J = ana_dev(mat,μ, T,rho,y,A,Z,m)
        J⁻¹ = pinv(J)
        #deter = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        #J⁻¹ = 1.0/deter * [J[2,2] -J[2,1]; -J[1,2] J[1,1]]
        #println(">>> fun", fun(μ))
        #println(">>> J⁻¹ ", J⁻¹)
        #println(det(J))
        μⁱ⁺¹ = μ .- (zaehler/50.0).*[J⁻¹[1,1]*fun(μ)[1] + J⁻¹[1,2]*fun(μ)[2]; J⁻¹[2,1]*fun(μ)[1] + J⁻¹[2,2]*fun(μ)[2]]
        μ = μⁱ⁺¹
        #println(">>> mu'", μ)
        #println(sqrt(f(F,μ,T,y,rho,A,Z,m)'f(F,μ,T,y,rho,A,Z,m)))
        #print("-----",μⁱ⁺¹,"\n")
        #ϵ = sum(mass_fraction([μⁱ⁺¹[1],μⁱ⁺¹[2]],T,rho,A,Z,m)) - 1.0
        ϵ = fun(μ)[1]^2 + fun(μ)[2]^2
        println(zaehler, "  ", ">>> ϵ >>>", ϵ," ", f(F,x,T,y,rho,A,Z,m))
    end
    println("iterations: ", zaehler)
    return μ
end




end
