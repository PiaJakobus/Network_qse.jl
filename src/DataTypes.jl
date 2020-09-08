using Base


"""
Mass [MeV]
--> M * [meverg] / c^2 = [erg]
"""
struct AtomicProperties
    #TODO: convert all to MeV
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
        #(A * m_u * c^2 / meverg + Δ) - ((Z * (m_p + m_e) + (A - Z) * m_n) * c^2 / meverg),
        Z, A, s, Δ, ω)
end

Base.show(io::IO, a::AtomicProperties) = print(io, "$(a.name), Z$(a.Z)A$(a.A), M = $(a.M), Eb = $(a.Eb), Δ = $(a.Δ)")
