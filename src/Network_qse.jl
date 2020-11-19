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



function testing(yrange::Vector, trange::Vector, rrange::Vector, Φ = 170.0)
    srange = Array{Float64, 1}(undef, size(trange, 1))
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.initial_guess(a)
    c12 = find_el("C12", a)
    for (i,y) in enumerate(yrange)
        for (j, t) in enumerate(trange)
            tmp = initial_guess(a)
            for (k, r) in enumerate(rrange)
                any(isnan.(tmp)) ? tmp = Network_qse.initial_guess(a) : nothing
                #r = Φ * (t/1.0e9)^3 * 1.0e5
                tmp = Network_qse.MultiNewtonRaphson(tmp, t, r, y, a)
                res[:,i,j,k] = Network_qse.x_i(tmp, t, r, a)
                srange[j] = sum(res[:,i,j,k][c12:end])
                println(">>>> ", j, " ", " sum ",sum(res[:,i,j,k]))
            end
        end
    end
    save("./NSE_table.jld", "data", res)
    save("./NSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", srange)
    open("./README_NSE.txt"; write=true) do f
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range, T-range, rho-range,  clust-Range\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", srange])
    end
    return res, srange
end


#TODO: Phi and percentage above NSE in MultiQSENR als Funktionsparameter uebergeben.
#TODO: percentage above NSE in Readme schreiben
function testing_QSE(yrange::Vector, trange::Vector, rrange::Vector, srange::Vector)
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 5}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1), size(srange, 1))
    tmp = Network_qse.qse_initial_guess(a)
    ig = tmp
    #q = 0.5
    Φ = 170
    for (i,y) in enumerate(yrange)
        for (j, t) in enumerate(trange)
            tmp = Network_qse.qse_initial_guess(a)
            for (k, r) in enumerate(rrange), (l, q) in enumerate(srange)
                any(isnan.(tmp)) ? tmp = Network_qse.qse_initial_guess(a) : nothing
                #r = Φ * (t/1.0e9)^3 * 1.0e5
                tmp = Network_qse.QSE_MultiNewtonRaphson(tmp, t, r, y, q, a)
                x_nse, x_qse = Network_qse.x_i_QSE(tmp, t, r, a)
                res[:,i,j,k,l] = vcat(x_nse,x_qse)
                #srange[j] = sum(x_qse)
                println(">>>> jl: ",j,"  ", l, " ", " sum X: ",sum(res[:,i,j,k,l]), "   sum X_cl: ", sum(x_qse))
                println("T    ", t, " ", " rrange    : ",r, "   sum X_cl: ", srange[l])
            end
        end
    end
    save("./QSE_table.jld", "data", res)
    save("./QSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", srange)
    open("./README_QSE.txt"; write=true) do f
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range, T-range, rho-range,  clust-Range\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", srange])
    end
    return res, srange
end



end
