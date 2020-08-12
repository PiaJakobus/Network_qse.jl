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
function initial_partition_function(ω,A,Z,s,m)::Array{Float64,2}
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

function mass_fraction(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    E_b =  (m .- Z*m_p .+ N*m_n).*const_meverg
    β = 1.0/(T*const_k_B)
    prefact = [exp.(inter_pr1[el](T)) for el in 1:length(A)]
    result = ((prefact / ρ) .* exp.((μ[2] .* Z .+ μ[1] .* N .- E_b).*β))
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
    β = 1.0/(const_k_B*T)
    #J[1,2] = sum(β.*Z.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
    #J[2,1] = sum((β.*N.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y))
    #J[2,2] = sum((β.*Z.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y))
    #J[1,1] = sum(β.*N.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
    J[1,1] = sum(β.*N.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
    J[1,2] = sum(β.*Z.*mass_fraction(μ, T,rho,A,Z,m))/sum(mass_fraction(μ, T,rho,A,Z,m))
    J[2,1] = sum((β.*N.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y))
    J[2,2] = sum((β.*Z.*Z./A).*mass_fraction(μ, T,rho,A,Z,m))/(sum((Z./A).*mass_fraction(μ, T,rho,A,Z,m))*log(y))
    #println("J[1,1]: ", J[1,1], " J[1,2]: ", J[1,2])
    #println("J[2,1]: ", J[2,1], " J[2,2]: ", J[2,2])
    #println(".... ",sum(mass_fraction(μ, T,rho,A,Z,m)))
    return J
end

function ana_dev!(M,μ, T,rho,y,A,Z,m)
    N = A .- Z
    β = 1.0/(const_k_B*T)
    M[1,1] = sum(β.*N.*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/sum(exp.(log_mass_fraction(μ, T,rho,A,Z,m)))
    M[1,2] = sum(β.*Z.*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/sum(exp.(log_mass_fraction(μ, T,rho,A,Z,m)))
    M[2,1] = sum((β.*N.*Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/(sum((Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))*log(y))
    M[2,2] = sum((β.*Z.*Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))/(sum((Z./A).*exp.(log_mass_fraction(μ, T,rho,A,Z,m)))*log(y))
end


J = ana_dev!(Array{Float64,2}(undef, (2,2)),[-1e-4, 1e-3],8e9,0.6,1e7,A,Z,m)
deter = J[1,1]*J[2,2] - J[1,2]*J[2,1]
det(J)
J⁻¹ = 1.0/deter * [J[2,2] -J[1,2]; -J[2,1] J[1,1]]
pinv(J)


"""return index in (A/Z/m/s)"""
function find_nucl(aa,zz,A,Z)
    y = findall(x->x == zz, Z)
    return y[findall(x->x==aa, A[y])]
end


tresh(min,tem,X) = map(j -> any(X[tem,:,1,j] .> min) && [X[tem,:,1,j],j], 1:length(Z))
filter_out(fmin,min,t) = permutedims(hcat(filter(i->i!=false,fmin(min,t)[:,1])...))


function my_newton_raphson(μ,T,rho,y,A,Z,m)
    J = zeros(Float64, 2,2)
    F = Array{Float64,2}(undef, 2, 1)
    fun(x) = f(F,x,T,y,rho,A,Z,m)
    N = A .- Z
    β = 1.0/(const_k_B*T)
    global ϵ = 1.0
    global zaehler = 0
    while ϵ > 1e-7
        zaehler += 1
        ana_dev!(J,μ, T,rho,y,A,Z,m)
        J⁻¹ = pinv(J)
        μⁱ⁺¹ = μ .- (zaehler/30.0).*[J⁻¹[1,1]*fun(μ)[1] + J⁻¹[1,2]*fun(μ)[2]; J⁻¹[2,1]*fun(μ)[1] + J⁻¹[2,2]*fun(μ)[2]]
        μ = μⁱ⁺¹
        ϵ = fun(μ)[1]^2 + fun(μ)[2]^2
        println(zaehler, "  ", ">>> ϵ >>>", ϵ)
    end
    println("iterations: ", zaehler)
    return μ
end






function solve(x_guess, f_root,dev_a,N_vec,i_solve)
    F             = Array{Float64,1}(undef,2)
    N, N_y, N_rho = N_vec
    sol_T         = Array{Float64,4}(undef,(N,N_y,N_rho,npart))
    chempot       = Array{Float64,4}(undef, (2,N,N_y,N_rho))
    sol_T_NR      = Array{Float64,4}(undef,(N,N_y,N_rho,npart))
    chempot_NR    = Array{Float64,4}(undef, (2,N,N_y,N_rho))
    t,rho,y       = 3.5e9,1e9,0.5
    i,j,k         = ones(Int64,3)
    #(i_solve == 0) && (solver(F,x,t,y,rho) = nlsolve((F,x)->f_root(F,x,t,y,rho,A,Z,m),x_guess,autodiff = :forward).zero)
    (i_solve == 1) && (solver(F,x,t,y,rho) = nlsolve((F,x)->f(F,x,t,y,rho,A,Z,m), (J,x)->ana_dev!(J,x,t,rho,y,A,Z,m), x_guess).zero)
    #(i_solve == 2) && (solver(F,x,t,y,rho) = my_newton_raphson(F,x,t,y,rho,A,Z,m))
    #println(i_solve)
    for (j,y) in enumerate(LinRange(0.42,0.5,N_y))
    #for (i,t) in enumerate(LinRange(2.5e8,12e9,N))
    #for (k,rho) in enumerate(LinRange(1e9,1e9,N_rho))
        sol = solver(F,x,t,y,rho)
        sol_NR = my_newton_raphson(F,x_guess,t,y,rho,A,Z,m)
        println(">>>>>>>>>>>>>>>>>>>")
        println("μₙ,μₚ:   ", sol)
        println("NR μₙ,μₚ:   ", sol_NR)
        println("T,y,rho: ", t, " ",y," ", rho)
        println("res_x:   ", sum(exp.(log_mass_fraction(sol,t,rho,A,Z,m)))-1)
        println("res_y:   ", logsumexp(log_charge_neutrality(sol,t,rho,A,Z,m))/log(y) -1)
        println("res_x_NR:   ", sum(exp.(log_mass_fraction(sol_NR,t,rho,A,Z,m)))-1)
        println("res_y_NR:   ", logsumexp(log_charge_neutrality(sol_NR,t,rho,A,Z,m))/log(y) -1)
        sol_T[i,j,k,:] = exp.(log_mass_fraction(sol, t,rho,A,Z,m))
        sol_T_NR[i,j,k,:] = exp.(log_mass_fraction(sol_NR, t,rho,A,Z,m))
        chempot[:,i,j,k] = sol
        chempot_NR[:,i,j,k] = sol_NR
        x_guess = sol
    end
    #end
    #end
    return sol_T,chempot,sol_T_NR,chempot_NR
end

x_guess = [1.3954920692562653e-5, -1.648082101646189e-5]
X_i,mu,x_NR,mu_NR = solve(x_guess,f,ana_dev,[1,60,2],1)





"""
-------------------------------- plotting ----------------------------
"""
plr = plot(LinRange(0.42,0.5,60),filter_out(tresh,1e-3,1)[1,1],yaxis=:log,label = "Z"*string(Int64(Z[filter_out(tresh,1e-3,1)[1,2]]))*", A"*string(Int64(A[filter_out(tresh,1e-3,1)[1,2]])))
specs = filter_out(tresh,0.001,1)
for (i,el) in enumerate(specs[:,1])
    plr = plot!(LinRange(0.42,0.5,60),el,label = "Z"*string(Int64(Z[specs[i,2]]))*", A"*string(Int64(A[specs[i,2]])), yaxis=:log,ylims=(1e-5,1))
end
plr


dict = Dict("n" => [1,0], "p" => [1,1], "he3" => [3,2], "o16" => [16,8], "ne20" => [20,10],
"si28" => [28,14], "ti50" => [50,22],"chr52" => [52,24], "fe54" => [54,26], "fe58" => [58,26],"cob55" => [55,27], "cop55" => [55,29],"ni56" => [56,28], "fe56" => [56,26],
"ni62" => [62,28])


range_T = LinRange(2.5e9,1.2e10,40)
vary_T = plot(range_T,X_i[:,1,1,find_nucl(1,0,A,Z)[1]] .+ X_i[:,1,1,find_nucl(1,1,A,Z)[1]] .+ 0.0000001,yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = "p+n")#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
for k in keys(dict_coco)
    plot!(range_T,X_i[:,1,1,find_nucl(dict_coco[k][1],dict_coco[k][2],A,Z)[1]] .+ 0.0000001,ylims =(1e-5), yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)
end
vary_T


range_y = LinRange(0.41,0.5,60)
vary_y = plot(range_y,X_i[1,:,1,find_nucl(52,26,A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "Yₑ", ylabel = "Xᵢ",label = "fe52")#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
for k in keys(dict)
    plot!(range_y,X_i[1,:,1,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001,yaxis=:log, legend=:right,xlabel = "Yₑ", ylabel = "Xᵢ",label = k)
end
vary_y
savefig(fig_auto, "vary_ye_all_rho1e7_t9e9.pdf")
savefig(vary_y,"vary_ye_t3.5e9_rho9e9.pdf")


range_rho = LinRange(1e7,1e10,N_rho)
vary_rho = plot(range_rho,sol_T_auto[:,1,1,find_nucl(1,0,A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = "p")#,yticks = ([1e-2,0.1,1], ["1e-2", "0.1","1"]))
for k in keys(dict)
    plot!(range_rho,X_i[1,1,:,find_nucl(dict[k][1],dict[k][2],A,Z)[1]] .+ 0.0000001, yaxis=:log,xlabel = "T [K]", ylabel = "Xᵢ",label = k)
end
vary_rho





end
