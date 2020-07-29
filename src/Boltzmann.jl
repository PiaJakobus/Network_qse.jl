module Boltzmann

#include("Io.jl")
include("Constants.jl")
export initial_partition_function

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
    ω,A,Z,s,m =  extract_partition_function()[1:5]
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


end
