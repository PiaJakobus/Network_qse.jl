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
    srange = Array{Float64, 1}(undef, size(trange, 1))
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 4}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1))
    tmp = Network_qse.initial_guess(a)
    c12 = find_el("C12", a)
    sp = Network_qse.StepParameter(-50,50,30)
    for (i,y) in enumerate(yrange)
        for (j, t) in enumerate(trange)
            tmp = initial_guess(a)
            for (k, r) in enumerate(rrange)
                th = Network_qse.ThermoProperties(t, r, y, -5.0)
                ff = Network_qse.Func(2, x -> Network_qse.nse_condition(x, th, a), x -> Network_qse.df_nse_condition(x,th,a), false)
                any(isnan.(tmp)) ? tmp = Network_qse.initial_guess(a) : nothing
                #tmp = Network_qse.MultiNewtonRaphson(tmp, th, a)
                tmp = Network_qse.MultiNewtonRaphson(tmp, ff, th, a, sp)
                res[:,i,j,k] = Network_qse.x_i(tmp, th, a)
                srange[j] = sum(res[:,i,j,k][c12:end])
                println(">>>> ", j, " ", " sum ",sum(res[:,i,j,k]))
            end
        end
    end
    save("./NSE_table.jld", "data", res)
    save("./NSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", srange)
    open("./README_NSE.txt"; write=true) do f
        write(f, "# Netwon Raphson parametrization: $(sp)\n")
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range, T-range, rho-range,  clust-Range\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", srange])
    end
    return res, srange
end



function testing_QSE(yrange::Vector, trange::Vector, rrange::Vector, srange::Vector)
    a = Network_qse.extract_partition_function()
    res = Array{Float64, 5}(undef, size(a,1), size(yrange, 1), size(trange, 1), size(rrange, 1), size(srange, 1))
    sp = Network_qse.StepParameter(-10,10,50)
    for (i,y) in enumerate(yrange)
        for (j, t) in enumerate(trange)
            tmp = Network_qse.qse_initial_guess(a)
            for (k, r) in enumerate(rrange), (l, q) in enumerate(srange)
                th = Network_qse.ThermoProperties(t, r, y, q)
                ff = Network_qse.Func(3, x -> Network_qse.qse_condition(x, th, a), x -> Network_qse.df_qse_condition(x,th,a), false)
                any(isnan.(tmp)) ? tmp = Network_qse.qse_initial_guess(a) : nothing
                tmp = Network_qse.MultiNewtonRaphson(tmp, ff, th, a, sp)
                #tmp = Network_qse.QSE_MultiNewtonRaphson(tmp, th, a)
                x_nse, x_qse = Network_qse.x_i_QSE(tmp, th, a)
                res[:,i,j,k,l] = vcat(x_nse,x_qse)
                #srange[j] = sum(x_qse)
                println(">>>> jl: ",j,"  ", l, " ", " sum X: ",sum(res[:,i,j,k,l]), "   sum X_cl: ", sum(x_qse))
                println("T    ", t, " ", " rrange    : ",r, "   sum X_cl: ", 1.0 - 10.0^(srange[l]))
            end
        end
    end
    save("./QSE_table.jld", "data", res)
    save("./QSE_params.jld", "trange", trange, "yrange", yrange, "rrange", rrange, "srange", 1.0 .- 10.0.^(srange))
    open("./README_QSE.txt"; write=true) do f
        write(f, "# Netwon Raphson parametrization: $(sp)\n")
        write(f, "# cgs units - rho = const = $(rrange[1]) g/cm3\n")
        write(f, "# load data in variable with load(\"data.jld\")[\"data\"]\n")
        write(f, "# load parameters with load(\"data.jld\")[\"yrange\"]\n")
        write(f, "# y-range, T-range, rho-range,  clust-Range\n")
        writedlm(f, [yrange, "\n",  trange, "\n", rrange, "\n", 1.0 .- 10.0.^(srange)])
    end
    return res, 1.0 .- 10.0.^(srange)
end



end
