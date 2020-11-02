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
        #println(ϵ)
    end
    #println(sum(x_i(x, T, rho, ap)), " QSE Cluster:  ", sum(x_i(x, T, rho, ap)[find_el("C12", ap):end]))
    return x
end


"""
    QSE_MultiNewtonRaphson(x::Vector, T, rho, y, ap)
∂/∂xᵢ [f₁,...,fₙ]
"""
function QSE_MultiNewtonRaphson(x::Vector, T, rho, y, x_qse, ap, ind = find_el("C12", ap))
    zaehler = 0
    ϵ = 1.0
    mu_nse = x[1:2]
    while abs(ϵ) > 1e-10
        df = Network_qse.df_qse_condition(x, T, rho, ap)
        mu_nse = MultiNewtonRaphson(mu_nse, T, rho, y, ap)
        x_qse = 0.8 * sum(x_i(mu_nse, T, rho, ap)[ind:end])
        #df1 = ForwardDiff.jacobian(x -> qse_condition(x, T, rho, y, x_qse, ap), x)
        f  = qse_condition(x, T, rho, y, x_qse, ap)
        #println("=========", x_qse)
        #detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
        #println(df .- df1)
        inv = inv_3x3(df)
        #inv = pinv(df)
        #println("--------------- ", inv)
        #println("----autodiff--- ", df1[1,:])
        #println(Network_qse.df_qse_condition(x, T, rho, ap))
        #println(x, T, rho)
        x = x .- min.(1, zaehler/50) * max.(min.(inv * f, 50), -50)
        ϵ = sqrt(f[1]^2 + f[2]^2 + f[3]^2)
        #TODO: keep the below line only for testing no need to call scr !!
        scr = Network_qse.x_i_QSE(x, T, rho, ap)
        println(zaehler, "  ", " >>> √ϵrror² >>> ", Float64(ϵ), "  sum x1 + x2:  ", sum(vcat(scr[1],scr[2])))
        zaehler += 1
    end
    return x
end
