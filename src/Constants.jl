const_m_B = 1.66e-24 # baryon mass
const_kmev = 8.61829e-11
const_meverg = 1.602e-6
const_ergmev = 1/const_meverg
const_k_B = 1.380658e-16
const_c = 2.99792458e10
const_h_barc = 197.327e-13
const_hh = const_h_barc / const_c * 2.0 * π * const_meverg


data_T = 10e9.*Float64[0.01, 0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]


#f(x) = saha_equation(x,1.1,2.2)
#h(x) = charge_neutrality(x,1.1,2.2)
##dxᵢ = x -> ForwardDiff.jacobian(f,x) # g = ∇f
#dyₑ = x -> ForwardDiff.jacobian(h,x)
#res1 = dxᵢ([1.1,2.2])
#res2 = dyₑ([1.1,2.2])


μ = [1.1,2.1]
##tt = saha_equation(μ,1.1,2.2)
#rr = sum.([saha_equation(μ, 1.9,1.1),saha_equation(μ, 1.9,1.1)])
#f_zeroth(μ) = sum.([saha_equation(μ, 1.9,1.1), charge_neutrality(μ, 1.9,1.1)]) - [1.0,2.0]
#f_zeroth(μ)
#optimize(f_zeroth, [0.0, 0.0],SimulatedAnnealing())
