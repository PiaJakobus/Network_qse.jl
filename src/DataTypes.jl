using Base

struct AtomicProperties
    name::String
    Z::Int64            # charge number
    A::Int64            # atomic number
    s::Float64          # spin
    m::Float64          # mass
    ω::Function         #\omega(T)
    AtomicProperties(Z::Int64, A::Int64, s::Float64, m::Float64, ω::Function) =
        new(Z==0 ? "1n" : PeriodicTable.elements[Z].symbol*string(A), Z, A, s, m, ω)
end

Base.show(io::IO, a::AtomicProperties) = print(io, "Atomic properties of $(a.name) with Z = $(a.Z)")
