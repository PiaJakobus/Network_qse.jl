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



function testing(yrange::Vector, trange::Vector, rrange::Vector, srange = Array{Float64, 1}(undef, 150))
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.initial_guess(a)
    c12 = find_el("C12", a)
    for (i,y) in enumerate(yrange), (j, t) in enumerate(trange), (k, r) in enumerate(rrange)
        tmp = Network_qse.MultiNewtonRaphson(tmp, t, r, y, a)
        res[:,i,j,k] = Network_qse.x_i(tmp, t, r, a)
        srange[j] = sum(res[:,i,j,k][c12:end])
        println(">>>> ", j, " ", " sum ",sum(res[:,i,j,k]))
    end
    save("./NSE_table.jld", "data", res)
    open("./README_NSE.txt"; write=true) do f
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# y-range, T-range\n")
        writedlm(f, [yrange, "\n",  trange])
    end
    return res, srange
end

function testing_QSE(yrange::Vector, trange::Vector, rrange::Vector, srange = Array{Float64, 1}(undef, 150))
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.qse_initial_guess(a)
    for (i,y) in enumerate(yrange), (j, t) in enumerate(trange), (k, r) in enumerate(rrange)
        tmp = Network_qse.QSE_MultiNewtonRaphson(tmp, t, r, y, 1.0, a)
        x_nse, x_qse = Network_qse.x_i_QSE(tmp, t, r, a)
        res[:,i,j,k] = vcat(x_nse,x_qse)
        srange[j] = sum(x_qse)
        println(">>>> ", j, " ", " sum ",sum(res[:,i,j,k]))
    end
    save("./QSE_table.jld", "data", res)
    open("./README_QSE.txt"; write=true) do f
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# y-range, T-range, clust-Range\n")
        writedlm(f, [yrange, "\n", trange, "\n", srange])
    end
    #close("./QSE_table.jld")
    #close("./README_QSE.txt")
    return res, srange
end



end
