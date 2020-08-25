"""
    prefactor(pf)

returns prefactor of X_i, as given in
http://cococubed.asu.edu/code_pages/nse.shtml
"""
function prefactor(pf)
    (T, rho) -> (pf.A * pf.ω(T) *
                (2 * pf.s + 1) *
                (2.0 * T * π * Network_qse.k_B * pf.M / Network_qse.hh^2)^1.5
                /(rho * Network_qse.N_A))
end


"""
    massCon(μ, T, rho, pf)

returns mass fraction for one element
"""
function massCon(μ::Vector, T::Float64, rho::Float64, pf)
    prefactor(pf)(T, rho) *
        exp((μ[2] * pf.Z + μ[1] * pf.A - pf.Eb) / (Network_qse.k_B*T))
end



"""
    chargeNeut(μ, T, rho, pf)

returns charge fraction for one element
"""
function chargeNeut(μ::Vector, T::Float64, rho::Float64, pf)
    prefactor(pf)(T, rho) * (pf.Z / pf.A) *
        exp((μ[2] * pf.Z + μ[1] * pf.A - pf.Eb) / (Network_qse.k_B*T))
end
