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
function NewtonRaphson(μ,T,rho,y,A,Z,m,pol)
    J = zeros(Float64, 2,2)
    F = Array{Float64,2}(undef, 2, 1)
    #fun(x) = f(x,T,y,rho,A,Z,m,pol)
    N = A .- Z
    β = 1.0/(const_kmev*T)
    global ϵ = 1.0
    global zaehler = 0
    while ϵ > 1e-10
        zaehler += 1
        ana_dev!(J,μ, T,rho,y,A,Z,m,pol)
        J⁻¹ = pinv(J)
        f!(F,μ,T,y,rho,A,Z,m,pol)
        μⁱ⁺¹ = μ .- [J⁻¹[1,1]*F[1] + J⁻¹[1,2]*F[2]; J⁻¹[2,1]*F[1] + J⁻¹[2,2]*F[2]]
        μ = μⁱ⁺¹
        ϵ = F[1]^2 + F[2]^2
        println(zaehler, "  ", ">>> ϵ >>>", ϵ,F)
    end
    println("iterations: ", zaehler)
    return μ
end
