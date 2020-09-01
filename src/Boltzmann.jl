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
    df_nse_condition!(J,μ, T,rho, ap)

computes Jacobian ∇f ∈ ℝ²×ℝ² with f ∈ ℝ², μ ∈ ℝ²
"""
function df_nse_condition!(dres,μ, T,ρ, y, ap)
    setprecision(4000) do
        res = zeros(size(μ))
        exp_i = ones(2)
        dexp_i = ones(BigFloat, 2,2)
        #dres = zeros(size(μ), size(μ))
        for apᵢ in ap
            pr_i = prefactor(apᵢ)(T, ρ)
            exp_i = exp(BigFloat((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev * T)))
            factor_i = [1.0, apᵢ.Z / apᵢ.A]
            res .= (pr_i .* factor_i .* exp_i) .+ res
            dexp_i = BigFloat.([apᵢ.Z apᵢ.A; apᵢ.Z*apᵢ.Z/apᵢ.A apᵢ.Z] ./ (k_B * T)) .* exp_i
            dres[1, 1] = BigFloat(dexp_i[1,1] / (pr_i .* factor_i[1] .* exp_i)) .+ dres[1, 1]
            dres[1, 2] = BigFloat(dexp_i[1,2] / (pr_i .* factor_i[1] .* exp_i)) .+ dres[1, 2]
            dres[2, 1] = BigFloat(dexp_i[2,1] / (pr_i .* factor_i[2] .* exp_i)) .+ dres[2, 1]
            dres[2, 2] = BigFloat(dres[1, 1])
        end
        dres
    end
end


"""
    nse_condition!(res, μ, T, ρ, y, ap; precision=400)
Mass conservation and charge neutrality
log (∑ᵢXᵢ) and log(∑ᵢ(Zᵢ/Aᵢ)Xᵢ / y)
"""
function nse_condition!(res, μ, T::Float64, ρ::Float64, y::Float64, ap; precision=4000)
    setprecision(precision) do
        #res = zeros(size(μ))
        for apᵢ in ap
            pr_i = prefactor(apᵢ)(T, ρ)
            exp_i = exp(BigFloat((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev*T)))
            factor_i = [1.0, apᵢ.Z / apᵢ.A]
            res .= (pr_i .* factor_i .* exp_i) .+ res
        end
        res[2] /= y
        res = log.(res)
    end
end



function exponent(μ::Vector, T::Float64, rho::Float64, apᵢ)
    (μ[2] * apᵢ.Z + μ[1] * (apᵢ.A - apᵢ.Z) - apᵢ.Eb) / (kmev*T)
end
