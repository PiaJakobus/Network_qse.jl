"""
    prefactor(pf)

returns prefactor of X_i, as given in
http://cococubed.asu.edu/code_pages/nse.shtml
"""
function prefactor(pf)
    (T, rho) -> (pf.A * pf.ω(T) *
                (2 * pf.s + 1) *
                (2.0 * T * π * Network_qse.k_B * pf.M * meverg / (c^2 * Network_qse.hh^2))^1.5
                /(rho * Network_qse.N_A))
end


"""
    massCon(μ, T, rho, ap)

returns mass fraction for one element
"""
function massCon(μ::Vector, T::Float64, rho::Float64, ap)
    prefactor(ap)(T, rho) *
        exp((μ[2] * ap.Z + μ[1] * ap.A - ap.Eb * Network_qse.meverg) / (Network_qse.k_B*T))
end


"""
    chargeNeut(μ, T, rho, ap)

returns charge fraction for one element
"""
function chargeNeut(μ::Vector, T::Float64, rho::Float64, ap)
    prefactor(ap)(T, rho) * (ap.Z / ap.A) *
        exp((μ[2] * ap.Z + μ[1] * ap.A - ap.Eb) / (Network_qse.k_B*T))
end


"""
    massCon_highPrec(μ, T, rho, ap)

returns mass fraction for one element with abritary precision
"""
function massCon_highPrec(μ1, μ2, T::Float64, ρ::Float64, ap)
    setprecision(40000) do
        BigFloat(prefactor(ap)(T, ρ)) *
        exp((BigFloat(μ1) * BigFloat(ap.Z) + BigFloat(μ2) * BigFloat(ap.A) - BigFloat(ap.Eb*meverg)) / BigFloat(k_B*T))
    end
end

"""
    massCon_highPrec(μ, T, rho, ap)

returns mass fraction for one element with abritary precision
"""
function chargeNeut_highPrec(μ1, μ2, T::Float64, ρ::Float64, ap)
    setprecision(40000) do
        BigFloat(prefactor(ap)(T, ρ)) * (BigFloat(ap.Z) / BigFloat(ap.A)) *
        exp((BigFloat(μ1) * BigFloat(ap.Z) + BigFloat(μ2) * BigFloat(ap.A) - BigFloat(ap.Eb*meverg)) / BigFloat(k_B*T))
    end
end
