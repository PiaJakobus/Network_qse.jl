"""
    logsumexp(arr)

max + log(sum(exp{arr - max}))
instead of
log(sum(exp(arr)))
"""
function logsumexp(arr::Vector)
    max = maximum(arr)
    dx = arr .- max
    sumexp = sum(exp.(dx))
    return max + log.(sumexp)
end

function find_el(el::String, ap::Any)
    filter(i -> (ap[i].name == el), 1:size(ap,1))[1]
end

function inv_3x3(m::Array)
    det = m[1,1] * (m[2,2] * m[3,3] - m[2,3] * m[3,2]) -
          m[1,2] * (m[2,1] * m[3,3] - m[2,3] * m[3,1]) +
          m[1,3] * (m[2,1] * m[3,2] - m[2,2] * m[3,1])
    m⁻¹ = (1.0/det) .*
        [m[2,2]*m[3,3]-m[2,3]*m[3,2] m[1,3]*m[3,2]-m[1,2]*m[3,3] m[1,2]*m[2,3]-m[1,3]*m[2,2];
         m[2,3]*m[3,1]-m[2,1]*m[3,3] m[1,1]*m[3,3]-m[1,3]*m[3,1] m[1,3]*m[2,1]-m[1,1]*m[2,3];
         m[2,1]*m[3,2]-m[2,2]*m[3,1] m[1,2]*m[3,1]-m[1,1]*m[3,2] m[1,1]*m[2,2]-m[1,2]*m[2,1]]
    return m⁻¹
end


function inv_2x2(df::Array)
    detInv = 1.0 / (df[1,1] * df[2,2] - df[1,2] * df[2,1])
    inv = ([df[2,2] -df[1,2]; -df[2,1] df[1,1]]) .* detInv
end

"""
    MultiNewtonRaphson(x::Vector, func::Any, th::Any, ap::Any, props::Any)
props contains parameters for NR. For NSE a min/max value of [-50,50] and
alpha = 30 is recommended. For QSE min/max = [-10,10] and alpha = 50
seems to run stable
"""
function MultiNewtonRaphson(x::Vector, func::Any, th::Any, ap::Any, props::Any)
    zaehler = 0
    ϵ = 1.0
    while abs(ϵ) > 1e-12
        f = func.f(x)
        inv = func.inv(x)
        x = x .- min.(1, zaehler/props.alpha) * max.(min.(inv * f, props.max), props.min)
        ϵ = sqrt(sum(f.^2))
        zaehler += 1
    end
    return x
end
