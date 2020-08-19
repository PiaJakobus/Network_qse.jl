include("src/Network_qse.jl")
Network_qse.theme(:juno)

ω,A,Z,s,m =  Network_qse.Io.extract_partition_function()

grids = [1,50,1]
reduced_netw = [[1,0],[1,1],[3,2],[56,28],[62,28]]
x_guess = [6.492724770496382e-6, -8.232974345747265e-6]
X_i,mu = Network_qse.solve(x_guess,Network_qse.f!,Network_qse.ana_dev!,grids,0)

X_i

function filt_x_i(g,min,X)
    if g[2]≠1
        println("y")
        tresh1(min,X) = map(j -> any(X[1,:,1,j] .> min) && [X[1,:,1,j],j], 1:length(X_i[1,1,1,:]))
        return permutedims(hcat(filter(i->i!=false,tresh1(min,X)[:,1])...))
    elseif g[1]≠1
        println("T")
        tresh2(min,X) = map(j -> any(X[:,1,1,j] .> min) && [X[:,1,1,j],j], 1:length(X_i[1,1,1,:]))
        return permutedims(hcat(filter(i->i!=false,tresh2(min,X)[:,1])...))
    elseif g[3]≠1
        println("rho!!")
        tresh3(min,X) = map(j -> any(X[1,1,:,j] .> min) && [X[1,1,:,j],j], 1:length(X_i[1,1,1,:]))
        return permutedims(hcat(filter(i->i!=false,tresh3(min,X)[:,1])...))
end
end




y_range = LinRange(0.42,0.5,50)
specs = filt_x_i(grids,1e-1,X_i)
paired = Network_qse.ColorBrewer.palette("Paired", length(specs[:,1]))



a = Network_qse.plot()
for (i,el) in enumerate(specs[:,1])
    x_pos = y_range[argmax(specs[i,1])] - 0.003
    y_pos = maximum(el)*(1 + 0.1)
    if Int(Z[specs[i,2]]) == 0
        period = "n"
    else
        period = Network_qse.elements[Int(Z[specs[i,2]])].symbol
    end
    neutr = string(Int(A[specs[i,2]]))
    a = Network_qse.plot!(y_range,el,lc=paired[i],#xaxis =:flip,
        yaxis=:log,ylims=(1e-3,1),legend=:false)
        Network_qse.annotate!(x_pos,y_pos ,period*neutr,Network_qse.Plots.font("Sans",9,paired[i]))
end
a
savefig(a,"replicated_cocub_y0_5rho1e7.pdf")
