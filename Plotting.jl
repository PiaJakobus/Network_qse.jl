include("src/Network_qse.jl")
Network_qse.theme(:juno)



ω,A,Z,s,m    =  Network_qse.Io.extract_partition_function()
grids        = [100,1,1]
reduced_netw = [[1,0],[1,1],[3,2],[4,2],[56,26],[56,28]]

x_new    = [-4.276576147076027e7, 6.0494486950703084e7]
x_guess0     = [-1.466278915407641e-4, -1.779705485617325e-4]
x_guess1     = [1.458606518719055e-5, -1.76605466448629e-5]
x_guess2     =  [6.492724770496382e-6, -8.232974345747265e-6]
t,rho,y      =  9e7,1e7,0.5
X_i,mu,vals  = Network_qse.solve(x_guess0,Network_qse.f!,Network_qse.ana_dev!,grids,0,t,rho,y,reduced_netw,false)

#[7.8150285728978e-6, -9.641632725249186e-6]

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


function plot_network(specs,y_range,X_i,vals,grids,reduced_netw,A,Z,test_network)
    if test_network
        paired = ["red","orange","green","yellow","white"]
        A = map(x->x[1], reduced_netw)
        Z = map(x->x[2], reduced_netw)
    else
         paired = 1#Network_qse.ColorBrewer.palette("Paired", length(specs[:,1]))
    end
    plr = Network_qse.plot()
    for (i,el) in enumerate(specs[:,1])
        x_pos = y_range[argmax(specs[i,1])] * 1.006
        y_pos = maximum(el)*1.1
        if Int(Z[specs[i,2]]) == 0
            period = "n"
        else
            period = Network_qse.elements[Int(Z[specs[i,2]])].symbol
        end
        neutr = string(Int(A[specs[i,2]]))
        plr = Network_qse.plot!(y_range,el,#c=paired[i],#xaxis =:flip,
            yaxis=:log,ylims=(1e-2,1),legend=:false,title="T/rho/y "*string(vals),titlefontsize=8)
        Network_qse.annotate!(x_pos,y_pos ,period*neutr)#,Network_qse.Plots.font("Sans",10,paired[i]))
    end
    return plr
end

specs = filt_x_i(grids,1e-2,X_i)
p = plot_network(specs,LinRange(1e9,10e9,100),X_i,vals,grids,reduced_netw,A,Z,false)
Network_qse.savefig(p,"varyT.png")


X_i
