module Network_qse

#"""
# https://docs.julialang.org/en/v1/manual/unicode-input/
#"""
__precompile__(false)


include("dependencies.jl")
include("IO.jl")
include("Constants.jl")
include("DataTypes.jl")
include("Boltzmann.jl")
include("Tools.jl")



function testing(yrange::Vector, trange::Vector, rrange::Vector)
    #println(".....")
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.initial_guess(a[863])
    for (i,y) in enumerate(yrange), (j, t) in enumerate(trange), (k, r) in enumerate(rrange)
        tmp = Network_qse.MultiNewtonRaphson(tmp, t, r, y, a)
        res[:,i,j,k] = Network_qse.x_i(tmp, t, r, a)
        println(">>>> ", i, " ", " sum ",sum(res[:,i,j,k]))
    end
    save("./NSE_table.jld", "data", res)
    open("./README.txt"; write=true) do f
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# y-range, T-range\n")
        writedlm(f, [yrange trange])
    end
    return res
end





end
