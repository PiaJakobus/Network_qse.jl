"""
    logsumexp(arr)
Useful for large exponents.
Uses the trick:
max(arr) + log(sum(exp{arr - max(arr)}))
instead of log(sum(exp(arr)))
"""
function logsumexp(arr)
    max = maximum(arr)
    dx = arr .- max
    sumexp = sum(exp.(dx))
    return max + log.(sumexp)
end



"""
    MultiNewtonRaphson(x::Vector, T, rho, y, ap)

Returns solution for proton and neutron chemical potentials.
Jacobian computet with ForwardDiff but can be changed to analytic function
df_nse_condition.
"""
function MultiNewtonRaphson(x::Vector, T, rho, y, ap)
    zaehler = 0
    ϵ = 1.0
    while abs(ϵ) > 1e-10
        #df = Network_qse.df_nse_condition(x,T,rho,y,ap)
        df = ForwardDiff.jacobian(x -> nse_condition(x,T,rho,y,ap), x)
        f  = nse_condition(x, T, rho, y, ap)
        detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
        inv = ([-df[2,2] df[1,2]; df[2,1] -df[1,1]] * f) .* detInv
        x = x .+ min.(1, zaehler/30) * max.(min.(inv, 50), -50)
        # old Version:
        #det = df[1,1] * df[2,2] - df[1,2] * df[2,1]
        #inv1 = - (df[2,2] * f[1] - df[1,2]*f[2]) / det
        #inv2 =   (df[2,1] * f[1] - df[1,1]*f[2]) / det
        #x[1] = x[1] + min(1, zaehler/30) * max(min(inv1, 50), -50)
        #x[2] = x[2] + min(1, zaehler/30) * max(min(inv2, 50), -50)
        ϵ = sqrt(f[1]^2 + f[2]^2)
        #println(zaehler, "  ", ">>> √ϵrror² >>>", Float64(ϵ), ":   ", x[1], ">>>>", x[2], " xi: ", sum(x_i(x, T, rho, ap)))
        zaehler += 1
    end
    #println("----------------------------------------------------------------")
    #println("----------------------------------------------------------------")
    #println("iterations: ", zaehler, "  ϵ: ", ϵ, " μ: ", x,  " sum: ", " T: ")
    #        T, " rho: ", rho, " y: ", y, " sum: ", sum(x_i(x, T, rho, ap)))
    return x
end
