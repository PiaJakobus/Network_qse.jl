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




#TODO: test this!!
"""

∂/∂xᵢ [f₁,...,fₙ] good guess: [-4.394094904641707, -12.915712928215058]
"""
function MultiNewtonRaphson(guess::Vector, T, rho, y, ap)
    zaehler = 0
    x = guess
    ϵ = 1.0
    while abs(ϵ) > 1e-11
        df = df_nse_condition(x, T, rho, y, ap)
        f  = nse_condition(x, T, rho, y, ap)
        det = df[1,1] * df[2,2] - df[1,2] * df[2,1]
        J⁻¹ = (1.0/det) * [df[2,2] -df[2,1];-df[1,2] df[1,1]]
        #J⁻¹ = pinv(Float64.(df))
        #println("J⁻¹ ", Float64.(J⁻¹), " --df-- ",Float64.(df), x)

        x[1] = x[1] - min(1, zaehler/30) * max(min((J⁻¹[1,1] * f[1] + J⁻¹[1,2] * f[2]), 50), -50)
        x[2] = x[2] - min(1, zaehler/30) * max(min((J⁻¹[2,1] * f[1] + J⁻¹[2,2] * f[2]), 50), -50)
        ϵ = abs(sqrt(f[1]^2 + f[2]^2))
        println(zaehler, "  ", ">>> ϵrror² >>>", Float64(f[1]), ":   ", Float64(f[2]), "::::", Float64(ϵ))
        zaehler += 1
    end
    println("iterations: ", zaehler)
    return x
end
