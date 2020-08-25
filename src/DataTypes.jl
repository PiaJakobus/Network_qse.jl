using Base

struct AtomicProperties
    #TODO: \Delta , M and Eb could be contructed here, compiler doesnt like
    name::String
    Z::Int64            # charge number
    A::Int64            # atomic number
    s::Float64          # spin
    Δ::Float64          # mass excess
    M::Float64          # nuclear mass
    Eb::Float64         # Binding energy
    ω::Function         #\omega(T)
    AtomicProperties(Z::Int64, A::Int64, s::Float64, Δ::Float64,
        M::Float64, Eb::Float64, ω::Function) =
        new(Z==0 ? "1n" : PeriodicTable.elements[Z].symbol*string(A),
        Z, A, s, Δ, M, Eb, ω)
end

Base.show(io::IO, a::AtomicProperties) = print(io, "Atomic properties of $(a.name) with Z = $(a.Z)")
