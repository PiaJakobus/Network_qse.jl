"""
    prefactor(pf)

returns prefactor of X_i, as given in
http://cococubed.asu.edu/code_pages/nse.shtml
"""
function prefactor(pf)
    th -> (pf.A * pf.ω(th.T) *
                (2 * pf.s + 1) *
                (2.0 * th.T * π * Network_qse.k_B * pf.M * Network_qse.meverg /
                (Network_qse.c^2 * Network_qse.hh^2))^1.5
                /(th.rho * Network_qse.N_A))
end

"""
    initial_guess(a, ap_ni56 = (k -> a[find_el(k, a)])("Fe56"),
                T = 2e9, rho = 1e7)
index: 438 or 863 or 762
inverse of saha equation with only one species
and μₚ = μₙ. returns μ
"""
function initial_guess(a, ap_ni56 = (k -> a[find_el(k, a)])("Fe56");
                T = 2e9, rho = 1e7)
    scr2 = rho / (m_u * ap_ni56.A)
    λ3 = sqrt(2.0 * pi * k_B * T * ap_ni56.A * m_u / hh^2.0)^3.0
    mu = ones(2) .* (kmev * T * log(scr2 / λ3) - 510.0) / ap_ni56.A
    return mu
end

"""
    qse_initial_guess(a, ap_ni56 = (k -> a[find_el(k, a)])("Fe56"),
                T = 2e9, rho = 1e7)

"""
function qse_initial_guess(a, ap_ni56 = (k -> a[find_el(k, a)])("Fe56");
                T = 2e9, rho = 1e7)
    scr2 = rho / (m_u * ap_ni56.A)
    λ3 = sqrt(2.0 * pi * k_B * T * ap_ni56.A * m_u / hh^2.0)^3.0
    mu = ones(3) .* (kmev * T * log(scr2 / λ3) - 510.0) / ap_ni56.A
    mu[3] =  28.0 * mu[1] - 236.533
    return mu
end



"""
    df_nse_condition!(J,μ, th, ap)

computes Jacobian ∇f ∈ ℝ²×ℝ² with f ∈ ℝ², μ ∈ ℝ²
"""
function df_nse_condition(μ::Vector, th::Any, ap::Any)
        sum_exp  = 0.0
        sum_expY = 0.0
        df11, df12, df21, df22 = zeros(Float64, 4)
        dres = zeros(Float64, 2, 2)
        for apᵢ in ap
            pr_i = prefactor(apᵢ)(th)
            exp_i = exp((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb)
                    / (kmev * th.T))
            sum_exp = (pr_i * exp_i) + sum_exp
            sum_expY = (pr_i * exp_i * apᵢ.Z / apᵢ.A) + sum_expY
            df11 = pr_i * exp_i * apᵢ.Z / (kmev * th.T) + df11
            df12 = pr_i * exp_i * (apᵢ.A - apᵢ.Z) / (kmev * th.T) + df12
            df21 = pr_i * exp_i * apᵢ.Z * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + df21
            df22 = pr_i * exp_i * (apᵢ.A - apᵢ.Z) * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + df22
        end
        dres[1, 1] = df11 / sum_exp #log (∑ᵢXᵢ)
        dres[1, 2] = df12 / sum_exp
        dres[2, 1] = - dres[1,1] * (sum_expY / sum_exp) + df21 / sum_exp
        dres[2, 2] = - dres[1,2] * (sum_expY / sum_exp) + df22 / sum_exp
        return dres
end


"""
    df_qse_condition!(J,μ, th, ap)

computes Jacobian ∇f ∈ ℝ³×ℝ³ with f ∈ ℝ³, μ ∈ ℝ³
filter(i -> (a[i].name == "C12"), 1:size(a,1)) = 7
"""
function df_qse_condition(μ::Vector, th::Any, ap::Any, E_si = -236.533, i_c12 = find_el("C12", ap))
        sum_exp  = 0.0
        sum_expY = 0.0
        sum_exp_si = 0.0
        sum_expY_si = 0.0
        df11, df12, df21, df22, df13, df31, df33, df23, df32 = zeros(Float64, 9)
        tmp1, tmp2 = zeros(Float64, 2)
        dres = zeros(Float64, 3, 3)
        for apᵢ in ap[1:i_c12]
            pr_i = prefactor(apᵢ)(th)
            exp_i    = exp((μ[1] * apᵢ.Z + μ[2] * (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev * th.T))
            sum_exp    = (pr_i * exp_i) + sum_exp
            sum_expY   = (pr_i * exp_i * apᵢ.Z / apᵢ.A) + sum_expY

            df11 = pr_i * exp_i * apᵢ.Z / (kmev * th.T) + df11
            df12 = pr_i * exp_i * (apᵢ.A - apᵢ.Z) / (kmev * th.T) + df12

            df21 = pr_i * exp_i * apᵢ.Z * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + df21
            df22 = pr_i * exp_i * (apᵢ.A - apᵢ.Z) * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + df22
        end
        for apᵢ in ap[i_c12+1:end]
            pr_i = prefactor(apᵢ)(th)
            exp_i_si = exp((μ[1] * (apᵢ.Z - 14) + μ[2] * ((apᵢ.A -apᵢ.Z) -14) + μ[3] - (apᵢ.Eb .+ E_si)) / (kmev * th.T))
            sum_exp_si = (pr_i * exp_i_si) + sum_exp_si
            sum_expY_si   = (pr_i * exp_i_si * apᵢ.Z / apᵢ.A) + sum_expY_si
            df13 = pr_i * exp_i_si / (kmev * th.T) + df13
            df23 = pr_i * exp_i_si * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + df23

            tmp1 = pr_i * exp_i_si * (apᵢ.Z - 14) * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + tmp1
            tmp2 = pr_i * exp_i_si * ((apᵢ.A - apᵢ.Z) - 14) * (apᵢ.Z / apᵢ.A) / (kmev * th.T) + tmp2

            df31 = pr_i * exp_i_si * (apᵢ.Z .- 14) / (kmev * th.T) + df31
            df32 = pr_i * exp_i_si * ((apᵢ.A - apᵢ.Z) - 14) / (kmev * th.T) + df32
            df33 = pr_i * exp_i_si / (kmev * th.T) + df33
        end
        dres[1, 1] = (df11 + df31) / (sum_exp + sum_exp_si)
        dres[1, 2] = (df12 + df32) / (sum_exp + sum_exp_si)
        dres[1, 3] = df33 / (sum_exp + sum_exp_si)

        dres[2, 1] = - dres[1,1] * ((sum_expY + sum_expY_si) / (sum_exp + sum_exp_si)) + (df21 + tmp1) / (sum_exp + sum_exp_si)
        dres[2, 2] = - dres[1,2] * ((sum_expY + sum_expY_si) / (sum_exp + sum_exp_si)) + (df22 + tmp2) / (sum_exp + sum_exp_si)
        dres[2, 3] = - dres[1,3] * ((sum_expY + sum_expY_si) / (sum_exp + sum_exp_si)) + df23/(sum_exp + sum_exp_si)

        dres[3, 1] = df31 / sum_exp_si
        dres[3, 2] = df32 / sum_exp_si
        dres[3, 3] = df33 / sum_exp_si
        return dres
end

"""
    qse_condition(μ::Array{Float64,1}, th::Any,
            ap::Any, E_si = -236.533, i_c12 = find_el("C12", ap))
Mass conservation and charge neutrality
log (∑ᵢXᵢ),  log(∑ᵢ(Zᵢ/Aᵢ)Xᵢ / y), log (∑ⱼ Xⱼ) where
j loops over the heavy Cluster
"""
function qse_condition(μ::Vector, th::Any,
            ap::Any, E_si = -236.533, i_c12 = find_el("C12", ap))
        res = zeros(Real,2) # Real here for AtoDiff to work
        res_si = zeros(Real, 2)
        sol = zeros(Real, 3)
        for (in, apᵢ) in enumerate(ap[1:i_c12])
            pr_i = prefactor(apᵢ)(th)
            exp_i = exp((μ[1] * apᵢ.Z + μ[2] *
                    (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev * th.T))
            factor_i = [1.0, apᵢ.Z / apᵢ.A]
            res .= (pr_i .* factor_i .* exp_i) .+ res
        end
        for (in, apᵢ) in enumerate(ap[i_c12+1:end])
            pr_i = prefactor(apᵢ)(th)
            exp_i_si = exp((μ[1] * (apᵢ.Z - 14) + μ[2] * ((apᵢ.A -apᵢ.Z) -14) + μ[3] - (apᵢ.Eb .+ E_si)) / (kmev * th.T))
            factor_i = [1.0, apᵢ.Z / apᵢ.A]
            res_si .= (pr_i .* factor_i .* exp_i_si) .+ res_si
        end
        sol[2] = (res[2] + res_si[2]) / (res[1] + res_si[1]) - th.y
        sol[1] = log(res[1] + res_si[1])
        #println(res[1], res_si[1])
        sol[3] = log(res_si[1] / th.x_qse)
        return sol
end





"""
    nse_condition(μ::Vector, th::Any, ap::Any)
Mass conservation and charge neutrality
log (∑ᵢXᵢ) and log(∑ᵢ(Zᵢ/Aᵢ)Xᵢ / y)
"""
function nse_condition(μ::Vector, th::Any, ap::Any)
        res = zeros(Real,2) # Real here for AtoDiff to work
        for (in, apᵢ) in enumerate(ap)
            pr_i = prefactor(apᵢ)(th)
            exp_i = exp((μ[1] * apᵢ.Z + μ[2] *
                    (apᵢ.A -apᵢ.Z) - apᵢ.Eb) / (kmev * th.T))
            factor_i = [1.0, apᵢ.Z / apᵢ.A]
            res .= (pr_i .* factor_i .* exp_i) .+ res
        end
        res[2] = res[2] / res[1] - th.y
        res[1] = log(res[1])
        return res
end


"""
    x_i_QSE(μ::Vector, th::Any, a::Any;
                E_si = -236.533, ind = find_el("C12", a))
Returns abundansed for NSE and 1 QSE clusters. The boundaries are
defaulted as function arguments
"""
function x_i_QSE(μ::Vector, th::Any, a::Any;
                E_si = -236.533, ind = find_el("C12", a))
    nse = map(apᵢ -> Network_qse.prefactor(apᵢ)(th) * exp((μ[1] * apᵢ.Z +
            μ[2] * (apᵢ.A - apᵢ.Z) - apᵢ.Eb) / (Network_qse.kmev * th.T)), a[1:ind])
    qse = map(apᵢ -> Network_qse.prefactor(apᵢ)(th) * exp((μ[1] * (apᵢ.Z -14)
            + μ[2] * ((apᵢ.A - apᵢ.Z) -14) + μ[3] - (apᵢ.Eb .+ E_si))
            /(Network_qse.kmev * th.T)), a[ind+1:end])
    return nse, qse
end


"""
    x_i(μ::Vector, T::Float64, ρ::Float64, a)
Returns particle abundances for single NSE cluster
"""
function x_i(μ::Vector, th::Any, a::Any)
    map(apᵢ -> Network_qse.prefactor(apᵢ)(th) * exp((μ[1] * apᵢ.Z + μ[2] *
            (apᵢ.A - apᵢ.Z) - apᵢ.Eb) / (Network_qse.kmev * th.T)), a)
end
