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

function linear_interpolation(xₐᵣᵣ, yₐᵣᵣ, x)
    """
    forward interpolation
    xₐᵣᵣ = [xᵢ,xᵢ₊₁]
    yₐᵣᵣ = [yᵢ,yᵢ₊₁]
    """
    y⁺ = yₐᵣᵣ[1] + (x - xₐᵣᵣ[1])*(yₐᵣᵣ[2] - yₐᵣᵣ[1])/(xₐᵣᵣ[2] - xₐᵣᵣ[1])
    return y⁺
end
