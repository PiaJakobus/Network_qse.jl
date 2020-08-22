using Base

struct AtomicProperties
    name::String
    Z::Int64            # charge number
    A::Int64            # atomic number
    s::Float64          # spin
    m::Float64          # mass
    AtomicProperties(Z, A, s, m) = new(PeriodicTable.elements[Z].symbol*string(A), Z, A, s, m)
end

Base.show(io::IO, a::AtomicProperties) = print(io, "Atomic properties of $(a.name) with Z = $(a.Z)")

struct PartitionFunction
    ap::AtomicProperties
    ω::Float64      # partition function
    T::Float64      # temperature
end
