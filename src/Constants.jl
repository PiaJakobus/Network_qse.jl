

const_m_B = 1.66e-24 # baryon mass
const_kmev = 8.61829e-11
const_meverg = 1.602e-6
const_ergmev = 1/const_meverg
const_k_B = 1.380658e-16
const_c = 2.99792458e10
const_N_A = 6.02214076e23
const_h_barc = 197.327e-13
const_hh = (const_h_barc / const_c) * 2.0 * π * const_meverg
data_T = 1e9.*Float64[0.01, 0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]
range = (49.0,19.0,19.0)
m_n = 8.071
m_p = 7.289



dict2 = Dict("n" => [1,0], "p" => [1,1], "he3" => [3,2], "o16" => [16,8], "ne20" => [20,10], "mg24" => [24,12],
"si28" => [28,14], "s32" => [32,16], "ar36" => [36,18], "ca40" => [40,20], "ti44" => [44,22], "ti50" => [50,22],"cr48" => [48,24],
"chr52" => [52,24], "fe52" => [52,26], "fe54" => [54,26], "cob55" => [55,27], "cop55" => [55,29],"ni56" => [56,28], "fe56" => [56,26],
"v52" => [52,23], "ni62" => [62,28])



dict = Dict("n" => [1,0], "p" => [1,1], "he3" => [3,2], "o16" => [16,8], "ne20" => [20,10],
"si28" => [28,14], "ti50" => [50,22],"chr52" => [52,24], "fe54" => [54,26], "cob55" => [55,27], "cop55" => [55,29],"ni56" => [56,28], "fe56" => [56,26],
"ni62" => [62,28])

dict_coco = Dict("fe54" => [54,26],"he3" => [3,2],"fe56" => [56,26])


"""

fe   = findall(x->x == 26, Z)
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

function solve_alt(x_guess, f_root,dev_a)

    N, N_y, N_rho                     = [40,1,1]
    sol_T,sol_T_NR,sol_T_auto         = [Array{Float64,4}(undef,(N,N_y,N_rho,npart)) for i ∈ 1:3]
    chempot, chempot_auto, chempot_NR = [Array{Float64,4}(undef, (2,N,N_y,N_rho)) for i ∈ 1:3]
    t,rho,y = 9e9,1e7,0.41
    i,j,k = ones(Int64,3)
    #(i_solve == 1) && solver(F,x) = nlsolve((F,x)->f_root(F,x,t,y,rho,A,Z,m),x_guess,autodiff = :forward)
    #(i_solve == 2) && solver(F,x) = nlsolve((F,x)->f(F,x,t,y,rho,A,Z,m), (J,x)->ana_dev(J,x,t,rho,y,A,Z,m), x_guess)
    #(i_solve == 3) && solver(F,x) = my_newton_raphson(x,t,y,rho,A,Z,m)

    #for (j,y) in enumerate(LinRange(0.41,0.6,N_y))#, (k,rho) in enumerate(LinRange(1e9,1e9,N_rho))
    for (i,t) in enumerate(LinRange(1e9,12e10,N))
        #sol_auto = nlsolve((F,x)->f_root(F,x,t,y,rho,A,Z,m),x_guess,autodiff = :forward)#, method = :newton)#iterations = 1000)
        sol_NR   = my_newton_raphson(F,x,t,y,rho,A,Z,m)
        #sol = nlsolve((F,x)->f(F,x,t,y,rho,A,Z,m), (J,x)->ana_dev(J,x,t,rho,y,A,Z,m), [6.492724770496382e-6, -8.232974345747265e-6])
        #println("chempot_auto: ", sol_auto.zero)#, "chempot_NR: ", sol_NR)
        println("chempot_NR: ", sol_NR, f(F,sol_NR,t,y,rho,A,Z,m))
        #println("chempot_ana: ", sol.zero)
        #println("temp: ", t, " << >> ","y_e ",y," ", sum(exp.(log_mass_fraction([sol_auto.zero[1],sol_auto.zero[2]],t,rho,A,Z,m)))-1,"  ", logsumexp(log_charge_neutrality([sol_auto.zero[1],sol_auto.zero[2]],t,rho,A,Z,m))/log(y) -1)
        #sol_T[i,j,k,:]      = exp.(log_mass_fraction([sol.zero[1],sol.zero[2]], t,rho,A,Z,m))
        sol_T_NR[i,j,k,:]     = mass_fraction(sol_NR, t,rho,A,Z,m)
        #sol_T_auto[i,j,k,:] = exp.(log_mass_fraction([sol_auto.zero[1],sol_auto.zero[2]], t,rho,A,Z,m))
        #chempot[:,i,j,k] = sol.zero
        #chempot_auto[:,i,j,k] = sol_auto.zero
        chempot_NR[:,i,j,k] = sol_NR
        x_guess = sol_NR#sol_auto.zero

    end
    return sol_T_auto,chempot_auto
end


X_i2,mu2 = solve_alt(x_guess,f,ana_dev)
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

"""
