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
function MultiNewtonRaphson(x::Vector, fun::Function, dfun::Function, ϵ)#μ,T,rho,y,A,Z,m,pol)
    zaehler = 0
    while ϵ > 1e-10
        J⁻¹ = pinv(dfun(x))
        x[1] = x[1] - (J⁻¹[1,1] * fun(x)[1] + J⁻¹[1,2] * fun(x)[2])
        x[2] = x[2] - (J⁻¹[2,1] * fun(x)[1] + J⁻¹[2,2] * fun(x)[2])
        #println("new guess ", J⁻¹)
        println(zaehler, "  ", ">>> ϵrror² >>>", x, sqrt(x[1]^2 + x[2]^2))
        zaehler += 1
    end
    println("iterations: ", zaehler)
    return x
end
