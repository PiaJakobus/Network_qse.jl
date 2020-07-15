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
