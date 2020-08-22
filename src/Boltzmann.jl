"""
    prefactor
returns prefactor of X_i, as
http://cococubed.asu.edu/code_pages/nse.shtml
ω
"""
function prefactor(pf::PartitionFunction)
    fp₀ = pf.ω * (2 * pf.ap.s + 1)
    λ = sqrt(hh^2 / (2.0*pf.T*π*k_B*(pf.ap.A*m_B + pf.ap.m*meverg/c^2)))
    pf.ap.A*m_B * fp₀/λ^3
end
