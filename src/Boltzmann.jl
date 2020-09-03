"""
    prefactor(pf)

returns prefactor of X_i, as given in
http://cococubed.asu.edu/code_pages/nse.shtml
"""
function prefactor(pf)
    (T, rho) -> (pf.A * pf.ω(T) *
                (2 * pf.s + 1) *
                (2.0 * T * π * Network_qse.k_B * pf.M * Network_qse.meverg /
                (Network_qse.c^2 * Network_qse.hh^2))^1.5
                /(rho * Network_qse.N_A))
end



"""
    initial_guess(ap_ni56)

inverse of saha equation with only one species
and μₚ = μₙ. returns μ
"""
function initial_guess(T, rho, ap_ni56)
    scr2 = rho / (m_u * ap_ni56.A)
    λ3 = sqrt(2.0 * pi * k_B * T * ap_ni56.A * m_u / hh^2.0)^3.0
    mu = (kmev * T * log(scr2 / λ3) - 510.0) / ap_ni56.A
end


"""
    df_nse_condition!(J,μ, T,rho, ap)

computes Jacobian ∇f ∈ ℝ²×ℝ² with f ∈ ℝ², μ ∈ ℝ²
"""
function df_nse_condition(μ, T,ρ, y, ap)
    #setprecision(4000) do
        sum_exp = zeros(2)
        df = zeros(2, 2)
        dres = zeros(2, 2)
        for apᵢ in ap
            pr_i = prefactor(apᵢ)(T, ρ)
            exp_i = exp((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev * T)) .* [1.0, apᵢ.Z / apᵢ.A]
            sum_exp .= (pr_i .* exp_i) .+ sum_exp
            df[1, 1] = pr_i * exp_i[1] * apᵢ.Z / (k_B * T) + df[1,1]
            df[1, 2] = pr_i * exp_i[2] * apᵢ.A / (k_B * T) + df[1,2]
            df[2, 1] = pr_i * exp_i[1] * apᵢ.Z*apᵢ.Z/apᵢ.A / (k_B * T) + df[2,1]
            df[2, 2] = pr_i * exp_i[2] * (apᵢ.Z*(apᵢ.A-apᵢ.Z) / apᵢ.A) / (k_B * T) + df[2,2]
        end
        dres[1, 1] = df[1,1] / sum_exp[1]
        dres[1, 2] = df[1,2] / sum_exp[2]
        dres[2, 1] = df[2,1] / sum_exp[1]
        dres[2, 2] = df[2,2] / sum_exp[2]
        return dres
    #end
end




"""
    nse_condition!(res, μ, T, ρ, y, ap; precision=400)
Mass conservation and charge neutrality
log (∑ᵢXᵢ) and log(∑ᵢ(Zᵢ/Aᵢ)Xᵢ / y)
"""
function nse_condition!(μ, T::Float64, ρ::Float64, y::Float64, ap; precision=4000)
    #setprecision(precision) do
        res = zeros(2)
        for apᵢ in ap
            pr_i = prefactor(apᵢ)(T, ρ)
            exp_i = exp((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev*T))
            factor_i = [1.0, apᵢ.Z / apᵢ.A]
            res .= (pr_i .* factor_i .* exp_i) .+ res
            #println(scr[1])
        end
        res[2] /= (y * res[1])
        res = log.(res)
        return res
    #end
end




function exponent(μ::Vector, T::Float64, rho::Float64, apᵢ)
    (μ[2] * apᵢ.Z + μ[1] * (apᵢ.A - apᵢ.Z) - apᵢ.Eb) / (kmev*T)
end
