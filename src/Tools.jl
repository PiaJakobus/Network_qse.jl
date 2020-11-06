"""
    logsumexp(arr)

max + log(sum(exp{arr - max}))
instead of
log(sum(exp(arr)))
"""
function logsumexp(arr)
    max = maximum(arr)
    dx = arr .- max
    sumexp = sum(exp.(dx))
    return max + log.(sumexp)
end

function find_el(el::String, ap)
    filter(i -> (ap[i].name == el), 1:size(ap,1))[1]
end

function inv_3x3(m::Array)
    det = m[1,1] * (m[2,2] * m[3,3] - m[2,3] * m[3,2]) -
          m[1,2] * (m[2,1] * m[3,3] - m[2,3] * m[3,1]) +
          m[1,3] * (m[2,1] * m[3,2] - m[2,2] * m[3,1])
    m⁻¹ = (1.0/det) *
        [m[2,2]*m[3,3]-m[2,3]*m[3,2] m[1,3]*m[3,2]-m[1,2]*m[3,3] m[1,2]*m[2,3]-m[1,3]*m[2,2];
         m[2,3]*m[3,1]-m[2,1]*m[3,3] m[1,1]*m[3,3]-m[1,3]*m[3,1] m[1,3]*m[2,1]-m[1,1]*m[2,3];
         m[2,1]*m[3,2]-m[2,2]*m[3,1] m[1,2]*m[3,1]-m[1,1]*m[3,2] m[1,1]*m[2,2]-m[1,2]*m[2,1]]
    return m⁻¹
end


"""
    MultiNewtonRaphson(x::Vector, T, rho, y, ap)
[-7.0010491,-9.1]
∂/∂xᵢ [f₁,...,fₙ] good guess: [-4.394094904641707, -12.915712928215058]
"""
function MultiNewtonRaphson(x::Vector, T, rho, y, ap)
    zaehler = 0
    ϵ = 1.0
    while abs(ϵ) > 1e-4
        df = Network_qse.df_nse_condition(x, T, rho, ap)
        f  = nse_condition(x, T, rho, y, ap)
        detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
        inv = ([-df[2,2] df[1,2]; df[2,1] -df[1,1]] * f) .* detInv
        x = x .+ min.(1, zaehler/30) * max.(min.(inv, 50), -50)
        ϵ = sqrt(f[1]^2 + f[2]^2)
        zaehler += 1
    end
    #println(sum(x_i(x, T, rho, ap)), " QSE Cluster:  ", sum(x_i(x, T, rho, ap)[find_el("C12", ap):end]))
    return x
end


"""
    QSE_MultiNewtonRaphson(x::Vector, T, rho, y, ap)
∂/∂xᵢ [f₁,...,fₙ]
γ ≡ X_cl / X (mass fraction in heavy cluster (without units))
"""
function QSE_MultiNewtonRaphson(x::Vector, T, rho, y, R, ap, ind = find_el("C12", ap))
    zaehler = 0
    ϵ = 1.0
    mu_nse = x[1:2]
    while abs(ϵ) > 1e-8
        df = Network_qse.df_qse_condition(x, T, rho, ap)
        #mu_nse = MultiNewtonRaphson(mu_nse, T, rho, y, ap)
        #x_nse = x_i(mu_nse, T, rho, ap)
        #α_qse = 0.99#*sum(x_nse[ind+1:end])/sum(x_nse)
        x_qse = R
        scr = Network_qse.x_i_QSE(x, T, rho, ap)
        #println("NSE X     ", sum(x_nse), "    NSE X_cl: ", sum(x_nse[ind+1:end]),
        #"    QSE X_cl  ", sum(scr[2]), "    QSE X ", sum([(scr...)...]))
        #df1 = ForwardDiff.jacobian(x -> qse_condition(x, T, rho, y, x_qse, ap), x)
        f  = qse_condition(x, T, rho, y, x_qse, ap)
        inv = inv_3x3(df)
        x = x .- min.(1, zaehler/50) * max.(min.(inv * f, 10), -10)
        ϵ = sqrt(f[1]^2 + f[2]^2 + f[3]^2)
        #TODO: keep the below line only for testing no need to call scr !!
        #scr = Network_qse.x_i_QSE(x, T, rho, ap)
        #println(mu_nse, df, x, T, rho)
        #println(zaehler, "  ", " >>> √ϵrror² >>> ", Float64(ϵ), "  sum x1 + x2:  ", sum(vcat(scr[1],scr[2])))
        zaehler += 1
    end
    return x
end
