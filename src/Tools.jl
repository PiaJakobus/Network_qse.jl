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
function MultiNewtonRaphson(xIter::Vector, fun::Function)#μ,T,rho,y,A,Z,m,pol)
    #J = zeros(Float64, 2,2)
    #F = Array{Float64,2}(undef, 2, 1)
    #fun(x) = f(x,T,y,rho,A,Z,m,pol)
    #N = A .- Z
    #β = 1.0/(const_kmev*T)
    #global ϵ = 1.0
    #global zaehler = 0
    #f(x) = f(x, args)
    #df(x)
    df(x) = x-> Network_qse.ForwardDiff.jacobian(fun, xIter)
    while ϵ > 1e-10
        #zaehler += 1
        #ana_dev!(J,μ, T,rho,y,A,Z,m,pol)

        J⁻¹ = pinv(df(x))
        #f!(F,μ,T,y,rho,A,Z,m,pol)
        xIter .-= J⁻¹.*fun(x)
        #μ = μⁱ⁺¹
        ϵ = xIter - J⁻¹.*fun(x)
        println(zaehler, "  ", ">>> ϵ >>>", ϵ)
    end
    println("iterations: ", zaehler)
    return xIter
end
