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
index: 438
inverse of saha equation with only one species
and μₚ = μₙ. returns μ
"""
function initial_guess(T, rho, ap_ni56)
    scr2 = rho / (m_u * ap_ni56.A)
    λ3 = sqrt(2.0 * pi * k_B * T * ap_ni56.A * m_u / hh^2.0)^3.0
    mu = (kmev * T * log(scr2 / λ3) - 510.0) / ap_ni56.A
    return mu .* ones(2)
end


"""
    df_nse_condition!(J,μ, T,rho, ap)

computes Jacobian ∇f ∈ ℝ²×ℝ² with f ∈ ℝ², μ ∈ ℝ²
"""
function df_nse_condition(μ, T,rho, y, ap)
    setprecision(4000) do
        sum_exp = zeros(BigFloat, 2)
        df = zeros(BigFloat, 2, 2)
        dres = zeros(BigFloat, 2, 2)
        for apᵢ in ap
            #pr_i = prefactor(apᵢ)(T, ρ)
            exp_i = exp(BigFloat((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev * T))) .* BigFloat.([1.0, apᵢ.Z / apᵢ.A])
            sum_exp .= BigFloat.(prefactor(apᵢ)(T,rho) * exp_i) .+ sum_exp
            df[1, 1] = BigFloat(prefactor(apᵢ)(T,rho) * exp_i[1] * apᵢ.Z / (kmev * T)) + df[1,1]
            df[1, 2] = BigFloat(prefactor(apᵢ)(T,rho) * exp_i[1] * apᵢ.A / (kmev * T)) + df[1,2]
            df[2, 1] = BigFloat(df[1, 1] * (apᵢ.Z / apᵢ.A)) + df[2,1]
            df[2, 2] = BigFloat(df[1, 2] * (apᵢ.Z / apᵢ.A)) + df[2,2]
        end
        dres[1, 1] = df[1,1] / sum_exp[1]
        dres[1, 2] = df[1,2] / sum_exp[1]
        dres[2, 1] = - df[1,1] * (sum_exp[2] / sum_exp[1]) + df[2,1] / sum_exp[1]
        dres[2, 2] = - df[1,2] * (sum_exp[2] / sum_exp[1]) + df[2,2] / sum_exp[1]
        return dres
    end
end




"""
    nse_condition!(res, μ, T, ρ, y, ap; precision=400)
Mass conservation and charge neutrality
log (∑ᵢXᵢ) and log(∑ᵢ(Zᵢ/Aᵢ)Xᵢ / y)
"""
function nse_condition(μ, T::Float64, ρ::Float64, y::Float64, ap; precision=4000)
    setprecision(precision) do
        res = zeros(BigFloat,2)
        for apᵢ in ap
            pr_i = BigFloat(prefactor(apᵢ)(T, ρ))
            exp_i = exp(BigFloat((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev*T)))
            factor_i = BigFloat.([1.0, apᵢ.Z / apᵢ.A])
            res .= BigFloat.((pr_i .* factor_i .* exp_i) .+ res)
            #println(scr[1])
        end
        res[2] = res[2] / res[1] - BigFloat(y)
        #res .= res .- [1.0, y]
        res[1] = log(res[1])
        return res
    end
end




function exponent(μ, T::Float64, apᵢ)
    (μ[2] * apᵢ.Z + μ[1] * (apᵢ.A - apᵢ.Z) - apᵢ.Eb) / (kmev*T)
end
