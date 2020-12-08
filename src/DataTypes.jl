using Base


"""
    AtomicProperties
Contains mass, binding energy, charge and atomic Number
spin, mass excess and temperature dependend interpolated
partition function of an element
"""
struct AtomicProperties
    name::String
    M::Float64          # nuclear mass
    Eb::Float64         # Binding energy
    Z::Int64            # charge number
    A::Int64            # atomic number
    s::Float64          # spin
    Δ::Float64          # mass excess
    ω::Function         #\omega(T)
    AtomicProperties(Z::Int64, A::Int64, s::Float64, Δ::Float64,
        ω::Function) =
        new(Z==0 ? "1n" : PeriodicTable.elements[Z].symbol*string(A),
        A * m_u * c^2 / meverg + Δ,
        Δ - (Z * Δₚ + (A - Z) * Δₙ),
        Z, A, s, Δ, ω)
end

Base.show(io::IO, a::AtomicProperties) = print(io, "$(a.name), Z$(a.Z)A$(a.A), M = $(a.M), Eb = $(a.Eb), Δ = $(a.Δ)")



struct ThermoProperties
    T::Float64
    rho::Float64
    y::Float64
    x_qse::Float64
    ThermoProperties(T::Float64, rho::Float64, y::Float64, x_qse::Float64) =
    new(T, rho, y, 1.0 - 10^(x_qse))
end
Base.show(io::IO, th::ThermoProperties) = print(io, "T = $(th.T) [K], ρ = $(th.rho) [g/cm³], y = $(th.y), X₂ = $(th.x_qse)")


"""
    Func
Stores all information for "universal" Newton Raphson method (2x2 or 3x3)
If FD is on, ForwardDiff computes Jacobian
Example: ff = Network_qse.Func(3, x -> Network_qse.qse_condition(x, th, a), x -> Network_qse.df_qse_condition(x,th,a), false)
"""
# https://discourse.julialang.org/t/finding-the-number-of-arguments-of-an-anonymous-function/24394/20
struct Func
    dims::Int
    f::Function
    df::Function
    FD::Bool
    inv::Function
    Func(dims::Int, f::Function, df::Function, FD::Bool) =
    new(dims, f, FD ? df = x -> ForwardDiff.jacobian(f, x) : df, FD, if (dims == 2) x -> inv_2x2(df(x)) elseif (dims == 3) x -> inv_3x3(df(x)) else  x-> inv(df(x)) end)
end
Base.show(io::IO, ff::Func) = print(io, "$(ff.dims == 2 ? "Rootfinding f ∈ ℝ²" : "f ∈ ℝ³"). AutoDiff is $(ff.FD ? "on" : "off")")


"""
    StepParameter
Specifies min/max stepwidth and alpha parameter for Newton Raphson method
Example: sp = Network_qse.StepParameter(-10,10,30)
"""
struct StepParameter
    min::Float64
    max::Float64
    alpha::Float64
end
Base.show(io::IO, sp::StepParameter) = print(io, "min.(1, count/$(sp.alpha)) * max.(min.(inv * f, $(sp.max), $(sp.min))")
