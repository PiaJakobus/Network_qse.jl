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
[-7.0010491,-9.1]
∂/∂xᵢ [f₁,...,fₙ] good guess: [-4.394094904641707, -12.915712928215058]
"""
function MultiNewtonRaphson(guess::Vector, T, rho, y, ap)
    zaehler = 0
    x = guess
    ϵ = 1.0
    #T = 2e9
    #rho = 1e6
    while abs(ϵ) > 1e-10
        df = Network_qse.df_nse_condition(x,T,rho,y,ap)
        f  = nse_condition(x, T, rho, y, ap)
        det = df[1,1] * df[2,2] - df[1,2] * df[2,1]
        inv1 = - (df[2,2] * f[1] - df[1,2]*f[2]) / det
        inv2 =   (df[2,1] * f[1] - df[1,1]*f[2]) / det
        x[1] = x[1] + min(1, zaehler/30) * max(min(inv1, 50), -50)
        x[2] = x[2] + min(1, zaehler/30) * max(min(inv2, 50), -50)
        ϵ = abs(sqrt(f[1]^2 + f[2]^2))
        #println(zaehler, "  ", ">>> ϵrror² >>>", Float64(ϵ), ":   ", x[1], ">>>>", x[2], " xi: ", sum(x_i(x, T, rho, ap)))
        zaehler += 1
    end
    println("----------------------------------------------------------------------------------------------------------------")
    println("iterations: ", zaehler, "  ϵ: ", ϵ, " μ: ", x,  " sum: ", sum(x_i(x, T, rho, ap)))
    println("----------------------------------------------------------------------------------------------------------------")
    return x
end


function testing(T, rho,a)
    yrange = LinRange(0.5,0.42,30)
    res = Array{Float64, 2}(undef, 2, size(yrange,1))
    tmp = Network_qse.initial_guess(a[863])
    for (i,y) in enumerate(yrange)
        #println(tmp," ", t," ", rho)
        #println(Network_qse.MultiNewtonRaphson(tmp, t, rho, 0.5, a))
        res[:,i] = Network_qse.MultiNewtonRaphson(tmp, T, rho, y, a)
        tmp = res[:,i]
       end
    return res
end


# He4: 6
# ni56: 762
# Ni58: 764
# Fe54: 660
# Fe55: 661
# fe56: 662
# Fe58: 664
# Co55: 710
# Cr52: 563
# Cr54: 565
# Ti50: 470

# res = Network_qse.testing(rho,a)
#f(i) = Network_qse.x_i(res[:,i], trange[i], rho, a)
# all = map(i -> f(i), 1:size(trange,1))
# all = hcat(all...)
# plot(trange, all[762,:], yaxis=:log, label = "Ni56")
# plot(yrange, all[762,:], yaxis=:log, ylim = (1e-3,1), xflip=true, label = "Ni56")
