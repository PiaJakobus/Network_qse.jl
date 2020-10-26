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
    while abs(ϵ) > 1e-10
        df = Network_qse.df_nse_condition(x, T, rho, ap)
        f  = nse_condition(x, T, rho, y, ap)
        detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
        inv = ([-df[2,2] df[1,2]; df[2,1] -df[1,1]] * f) .* detInv
        x = x .+ min.(1, zaehler/30) * max.(min.(inv, 50), -50)
        ϵ = sqrt(f[1]^2 + f[2]^2)
        zaehler += 1
    end
    return x
end


"""
    QSE_MultiNewtonRaphson(x::Vector, T, rho, y, ap)
∂/∂xᵢ [f₁,...,fₙ]
"""
function QSE_MultiNewtonRaphson(x::Vector, T, rho, y, x_qse, ap)
    zaehler = 0
    ϵ = 1.0
    while abs(ϵ) > 1e-10
        #df = Network_qse.df_qse_condition(x, T, rho, ap)
        df = ForwardDiff.jacobian(x -> qse_condition(x, T, rho, y, x_qse, ap), x)
        f  = qse_condition(x, T, rho, y, x_qse, ap)
        #detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
        inv = inv_3x3(df)
        #println(Network_qse.df_qse_condition(x, T, rho, ap))
        #println(x, T, rho)
        x = x .+ min.(1, zaehler/10) * max.(min.(inv, 100), -100)
        ϵ = sqrt(f[1]^2 + f[2]^2)
        println(zaehler, "  ", ">>> √ϵrror² >>>", Float64(ϵ), ":   ", x[1], ">>>>", x[2], ">>>>", x[3], " xi: ", sum(x_i(x, T, rho, ap)))
        zaehler += 1
    end
    return x
end
