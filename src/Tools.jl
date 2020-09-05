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

∂/∂xᵢ [f₁,...,fₙ]
"""
function MultiNewtonRaphson(guess::Vector, fun::Function, dfun::Function, ϵ)#μ,T,rho,y,A,Z,m,pol)
    zaehler = 0
    x = guess
    while ϵ > 1e-10
        df = dfun(x)
        f = fun(x)
        J⁻¹ = pinv(df)
        x[1] = x[1] - (zaehler/40) * (J⁻¹[1,1] * f[1] + J⁻¹[1,2] * f[2])
        x[2] = x[2] - (zaehler/40) * (J⁻¹[2,1] * f[1] + J⁻¹[2,2] * f[2])
        #println("new guess ", J⁻¹)
        println(zaehler, "  ", ">>> ϵrror² >>>", f[1], ":   ", f[2], "   ", x)
        zaehler += 1
    end
    println("iterations: ", zaehler)
    return x
end
