module Network_qse

#"""
# https://docs.julialang.org/en/v1/manual/unicode-input/
#"""
__precompile__(false)


include("dependencies.jl")
include("IO.jl")
include("Constants.jl")
include("DataTypes.jl")
include("Boltzmann.jl")
include("Tools.jl")



function solve_var_T(N,rho, y, a)
    N_part = size(a,1)
    μ = initial_guess(a[762])
    chemPot = Array{Float64, 2}(undef, N,2)
    x_i = Array{Float64}(undef, N, N_part)
    for (i,t) in enumerate(LinRange(1e8,1e9,N))
        chemPot[i,:] = MultiNewtonRaphson(μ, t, rho, y, a)
        println("AA", chemPot[i,:])
        x_i[i,:] = map(apᵢ -> Network_qse.prefactor(apᵢ)(t, rho) * exp((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A - apᵢ.Z) - apᵢ.Eb) / (kmev*t)), a)
        println("BB")
        #map(a_i -> Network_qse.x_i(chemPot[i,:], t, rho, a_i), a)
        # sum(map(i -> Network_qse.x_i(μ,1e6,T,a[i]), 1:size(a,1)))
        println(">>>>>>>>>>>>>>>>>>> i: ", i, " T: ", t, " ∑ᵢ: ", chemPot[i,:], sum(x_i[i,:]))
    end
    return x_i, chemPot
end


end
