"""
    prefactor(pf)

returns prefactor of X_i, as given in
http://cococubed.asu.edu/code_pages/nse.shtml
"""
function prefactor(T, pf::AtomicProperties)
    #TODO: mit Stift und Papier testen
    fp₀ = pf.ω(T) * (2 * pf.s + 1)
    λ = sqrt(hh^2 / (2.0*T*π*k_B*(pf.A*m_B + pf.m*meverg/c^2)))
    pf.A*m_B * fp₀/λ^3
end
