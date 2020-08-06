

const_m_B = 1.66e-24 # baryon mass
const_kmev = 8.61829e-11
const_meverg = 1.602e-6
const_ergmev = 1/const_meverg
const_k_B = 1.380658e-16
const_c = 2.99792458e10
const_h_barc = 197.327e-13
const_hh = const_h_barc / const_c * 2.0 * π * const_meverg
data_T = 1e9.*Float64[0.01, 0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]
range = (49.0,19.0,19.0)
m_n = 8.071
m_p = 7.289


"""

fe   = findall(x->x == 26, Z)
chr  = findall(x->x == 24, G_all[3])
cob  = findall(x->x == 27, G_all[3])
ni   = findall(x->x == 28, G_all[3])
cop  = findall(x->x == 29, G_all[3])
ti   = findall(x->x == 22, G_all[3])

fe56  = fe[findall(x->x==56, A[fe])]
fe54  = fe[findall(x->x==50, A[findall(x->x == 26, G_all[3])])]
chr52 = chr[findall(x->x==52, A[findall(x->x == 24, G_all[3])])]
cob55 = cob[findall(x->x==55, A[findall(x->x == 27, G_all[3])])]
ni56  = ni[findall(x->x==56, A[findall(x->x == 28, G_all[3])])]
cop55 = cop[findall(x->x==55, A[findall(x->x == 29, G_all[3])])]
ti50  = ti[findall(x->x==50, A[findall(x->x == 22, G_all[3])])]



"""
