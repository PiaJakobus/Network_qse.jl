module Network_qse

#"""
# https://docs.julialang.org/en/v1/manual/unicode-input/
#"""
__precompile__(false)


#using Optim
using ForwardDiff
using Interpolations
using Plots
using NLsolve
#using Dierckx
using LinearAlgebra
#using IntervalRootFinding
using ColorBrewer
using ColorTypes
using PeriodicTable
using PlotThemes
include("Io.jl")
include("Constants.jl")


#export initial_partition_function
#export index_reduced_nse
#export nse_reduced
#export charge_neutrality
#export mass_conservation
#export f!
#export ana_dev!
export solve



"""
returns prefactor of X_i, as
http://cococubed.asu.edu/code_pages/nse.shtml
"""
function initial_partition_function(ω,A,Z,s,m)::Array{Float64,2}
    n_A = 1.0/const_m_B
    root_T⁻¹ = .√(1.0./(data_T))
    fp₀ = ω.*(2 .*s .+ 1)
    λ₀ = .√(const_hh^2/(2.0*π*const_k_B*(A*const_m_B .+ m*const_meverg/const_c^2)))
    λ = root_T⁻¹*λ₀
    prefac = vcat(map(i->((A[i]/n_A)*fp₀[i]/(λ[i].^3.0)/n_A), 1:length(fp₀))...)#
    return prefac
end





function index_reduced_nse(A,Z,A_i,Z_i)
    Z_arr = findall(x->x == Z_i, Z)
    iso = Z_arr[findall(x->x==A_i, A[findall(x->x == Z_i, Z)])]
    return iso
end


function nse_reduced(liste_el,arr,A,Z)
    return vcat(map(j -> arr[index_reduced_nse(A,Z,j[1],j[2])], liste_el)...)
end




function charge_neutrality(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector,pol)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    E_b =  (m .- Z*const_m_p .+ N*const_m_n).*const_meverg
    β = 1.0/(const_k_B*T)
    prefact = [exp.(pol[el](T)) for el in 1:length(A)]
    result = (prefact.*(Z./A)./ ρ).*exp.((μ[2] .* Z .+ μ[1] .* N .- E_b).*β)
    return result#(μ[2] .* Z .+ μ[1] .* N .- E_b).*β
end


function mass_conservation(μ::Vector,T::Float64,ρ::Float64,A::Vector,Z::Vector,m::Vector,pol)
    N = A .- Z
    result = zeros(eltype(μ),length(A))
    E_b =  (m .- Z*const_m_p .+ N*const_m_n)*const_meverg
    β = 1.0/(const_k_B*T)
    prefact = [exp.(pol[el](T)) for el in 1:length(A)]
    result = (prefact./ ρ).*exp.((μ[2] .* Z .+ μ[1] .* N .- E_b).*β)
    return result
end




function f!(F,x,T,yₑ,ρ,A,Z,m,pol)
    F[1] = sum(mass_conservation(x,T,ρ,A,Z,m,pol)) - 1
    F[2] = sum(charge_neutrality(x, T,ρ,A,Z,m,pol)) - yₑ
    #return F
end




function ana_dev!(J,μ, T,rho,y,A,Z,m,pol)
    N = A .- Z
    β = 1.0/(const_k_B*T)
    J[1,1] = sum(β.*N.*mass_conservation(μ, T,rho,A,Z,m,pol))
    J[1,2] = sum(β.*Z.*mass_conservation(μ, T,rho,A,Z,m,pol))
    J[2,1] = sum((β.*N.*Z./A).*mass_conservation(μ, T,rho,A,Z,m,pol))
    J[2,2] = sum((β.*Z.*Z./A).*mass_conservation(μ, T,rho,A,Z,m,pol))
end





function solve(x, f_root,dev_a,N_vec,i_solve,reduce = false)
    ω,A,Z,s,m =  Io.extract_partition_function()
    reduce && ((ω,A,Z,s,m) = map(x->nse_reduced(reduced_netw,x,A,Z),[ω,A,Z,s,m]))
    pr = initial_partition_function(ω,A,Z,s,m)
    linear_int = [LinearInterpolation(data_T, log.(pr[j,:])) for j in 1:length(A)]
    F             = Array{Float64,1}(undef,2)
    N, N_y, N_rho = N_vec
    sol_T         = Array{Float64,4}(undef,(N,N_y,N_rho,length(A)))
    chempot       = Array{Float64,4}(undef, (2,N,N_y,N_rho))
    sol_T_NR      = Array{Float64,4}(undef,(N,N_y,N_rho,length(A)))
    chempot_NR    = Array{Float64,4}(undef, (2,N,N_y,N_rho))
    t,rho,y       = 3.5e9,1e7,0.5
    i,j,k         = ones(Int64,3)
    (i_solve == 0) && (solver(F,x,t,y,rho) = nlsolve((F,x)->f_root(F,x,t,y,rho,A,Z,m,linear_int),x,autodiff = :forward).zero)
    #(i_solve == 1) && (solver(F,x,t,y,rho) = nlsolve((F,x)->f(F,x,t,y,rho,A,Z,m,linear_int), (J,x)->ana_dev!(J,x,t,rho,y,A,Z,m), x_guess).zero)
    #(i_solve == 2) && (solver(F,x,t,y,rho) = my_newton_raphson(x,t,y,rho,A,Z,m,linear_int))
    #(i_solve == 3) && (solver(x,t,y,rho) = IntervalRootFinding.roots(x->f(x,t,y,rho,A,Z,m,linear_int), x->ana_dev(x,t,rho,y,A,Z,m), X×Y, IntervalRootFinding.Newton))
    for (j,y) in enumerate(LinRange(0.42,0.5,N_y))
    #for (i,t) in enumerate(LinRange(2.5e7,10e9,N))
    #for (k,rho) in enumerate(LinRange(1e9,1e9,N_rho))
        #sol_NR = my_newton_raphson(F,sol,t,y,rho,A,Z,m)
        sol = solver(F,x,t,y,rho)
        println(">>>>>>>>>>>>>>>>>>>")
        println("μₙ,μₚ:   ", sol)
        #println("NR μₙ,μₚ:   ", sol_NR)
        println("T,y,rho: ", t, " ",y," ", rho)
        println("res_x:   ", sum(mass_conservation(sol,t,rho,A,Z,m,linear_int))-1)
        println("res_y:   ", sum(charge_neutrality(sol,t,rho,A,Z,m,linear_int)) - y)
        #println("res_x_NR:   ", sum(exp.(mass_conservation(sol_NR,t,rho,A,Z,m,linear_int)))-1)
        #println("res_y_NR:   ", logsumexp(charge_neutrality(sol_NR,t,rho,A,Z,m,linear_int))/log(y) -1)
        sol_T[i,j,k,:] = mass_conservation(sol, t,rho,A,Z,m,linear_int)
        #sol_T_NR[i,j,k,:] = exp.(mass_conservation(sol_NR, t,rho,A,Z,m,linear_int))
        chempot[:,i,j,k] = sol
        x = sol
        #chempot_NR[:,i,j,k] = sol_NR
    end
    #end
    #end
    return sol_T,chempot
end







end
