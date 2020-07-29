module TimeIt

if VERSION >= v"0.7.0"
    using Printf
end

export @timeit

macro timen(ex, n)
    quote
        local t0 = time_ns()
        for i = 1:$(esc(n))
            local val = $(esc(ex))
        end
        local t1 = time_ns()
        (t1 - t0) / 1.e9
    end
end

macro timeit(ex)
    quote
        local val = $(esc(ex))  # Warm up
        t = zeros(3)

        # Determine number of loops so that total time > 0.1s.
        n = 1
        for i = 0:9
            n = 10^i
            t[1] = @timen $(esc(ex)) n
            if t[1] > 0.1
                break
            end
        end

        # Two more production runs.
        for i = 2:3
            t[i] = @timen $(esc(ex)) n
        end
        best = minimum(t) / n

        # Format to nano-, micro- or milliseconds.
        factor = 1.
        if best < 1e-6
            factor = 1e9
            pre = "n"
        elseif best < 1e-3
            factor = 1e6
            pre = "\u00b5"
        elseif best < 1.
            factor = 1e3
            pre = "m"
        else
            factor = 1.
            pre = ""
        end
        @printf "%d loops, best of 3: %4.2f %ss per loop\n" n best*factor pre
        best
    end
end
a = zeros(Float64,10000,10000)
@timeit (a.^2).^2
@timeit a.^4
# loops, best of 3: 172.43 ms per loop
# loops, best of 3: 1.11 s per loop

function benchmark()
    a = zeros(Float64,10000,10000)
    benchmark_4(a::Array{Float64,2}) = a.^4
    benchmark_22(a::Array{Float64,2}) = (a.^2).^2
    benchmark_4(a)
    @timeit benchmark_4(a)
    benchmark_22(a)
    @timeit benchmark_22(a)
end

benchmark()
benchmark()

#  loops, best of 3: 170.89 ms per loop


end  # module
