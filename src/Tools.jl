function findnearest(a,x)
       n = length(a)
       n > 0 || return 0:-1
       i1 = searchsortedlast(a,x)
       i0 = i1
       if i1>0
           while i0>1 && a[i0-1]==a[i0]
               i0 -= 1
           end
           d = x-a[i1]
       else
           i0 = i1+1
           d = a[i1+1]-x
       end
       i2 = i1
       if i2<n && a[i2+1]-x<d
           i0 = i2+1
           d = a[i2+1]-x
           i2 += 1
       end
       while i2<n && a[i2+1]-x==d
           i2 += 1
       end
       return i0:i2
end

searchsortedlast([2,3.1,6],3)
# https://docs.julialang.org/en/v1/base/sort/



"""
    linear_interpolation(xₐᵣᵣ,yₐᵣᵣ,x)
find f(x) asuming f(x) = x, x∈[xᵢ,xᵢ₊₁]
xₐᵣᵣ = [xᵢ,xᵢ₊₁]
yₐᵣᵣ = [yᵢ,yᵢ₊₁]
"""
function linear_interpolation(xₐᵣᵣ::Array{Float64}, yₐᵣᵣ::Array{Float64}, x::Float64)
    y⁺ = yₐᵣᵣ[1] + (x - xₐᵣᵣ[1])*(yₐᵣᵣ[2] - yₐᵣᵣ[1])/(xₐᵣᵣ[2] - xₐᵣᵣ[1])
    return y⁺
end



"""
    Multivariate Newton raphson()
[xⁱ⁺¹₁..xⁱ⁺¹ₙ] = [xⁱ₁..xⁱₙ] - J⁻¹[f¹(xⁱ₁)..fⁿ(xⁱₙ)]
specificly:
calculate Hessian 2x2
[df[1]/dmun df[1]/dmup; df[2]/dumun df[2]/dmup]
J^-1 = 1/(ad-bc) * [d -b; -c a]
[mun,mup]' = [mun,mup] - J^-1 * f(mun,mup)
dXdμₙ   dXdμₚ
dYₑdμₙ  dYₑdμₚ
"""
function my_newton_raphson(μ::Array{Float64},T::Float64,ρ::Float64)
    y = 0.49
    G = extract_partition_function()
    f(μ::Array{Float64}) = sum.(saha_equation(μ, T, ρ)) - [1, y]
    #G = extract_partition_function()
    A, Z = G[2:3]
    β = const_k_B * T
    ϵ = 1.0
    while ϵ > 0.1
        println("error: ", ϵ)
        dXdμₙ  = sum((A .- Z)*β.*saha_equation(μ, T, ρ)[1])
        dXdμₚ  = sum(β*Z.*saha_equation(μ, T, ρ)[1])
        dYₑdμₙ = sum(dXdμₙ./(A.*Z))
        dYₑdμₚ = sum(dXdμₚ./(A.*Z))
        det = (dXdμₙ*dYₑdμₚ - dXdμₚ*dYₑdμₙ)
        J⁻¹ = 1.0/det * [dXdμₙ dXdμₚ; dYₑdμₙ dYₑdμₚ]
        μⁱ⁺¹ = μ .- J⁻¹.*f(μ)
        ϵ = √(abs(sum((μⁱ⁺¹ - μⁱ⁺¹).*(μⁱ⁺¹ - μⁱ⁺¹))))
    end
    return μ
end
sol = my_newton_raphson([1.1,2.1],2.1,2.2)
