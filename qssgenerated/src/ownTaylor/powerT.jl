# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#



#= The following method computes `a^float(n)` (except for cases like
Taylor0{Interval{T}}^n, where `power_by_squaring` is used), to
use internally `pow!`.
=#
powerT(a::Taylor0, n::Integer,cache1::Taylor0) = powerT(a, float(n),cache1)



#= 
#this is a power that want more caches and it is for taylor of ints which we do not use in qss...so maybe delete and use 1 cache
#the other two cases that want 2 caches are taylor^taylor and taylor^complex which we probably rarely encounter so make them use one cache
# in these rare cases it will be some allocs. it is ok.
    function powerT(a::Taylor0{T}, n::Integer,cache1::Taylor0,cache2::Taylor0) where {T<:Integer}
        n == 0 && return one(a,cache1)
        n == 1 && return a#copy(a)
        n == 2 && return square(a,cache1)
        n < 0 && throw(DomainError())
        return power_by_squaring(a, n,cache1,cache2)
    end

   function powerT(a::Taylor0{Rational{T}}, n::Integer,cache1::Taylor0,cache2::Taylor0) where {T<:Integer}
        n == 0 && return one(a,cache1)
        n == 1 && return a
        n == 2 && return square(a,cache1)
        n < 0 && return inv( powerT(a,-n,cache1,cache2) )
        return power_by_squaring(a, n,cache1,cache2)
    end

    powerT(a::Taylor0, x::Rational,cache1::Taylor0,cache2::Taylor0) = powerT(a,x.num/x.den,cache1,cache2)

     function powerT(a::Taylor0, b::Taylor0,cache1::Taylor0,cache2::Taylor0) 
         *( b,log(a,cache1),cache2 ) #this returns cache2
         cache1.coeffs.=0.0
         exp(cache2,cache1)
         
     end

    function powerT(a::Taylor0, x::T,cache1::Taylor0,cache2::Taylor0) where {T<:Complex} 
        *( x,log(a,cache1),cache2 ) #this returns cache2
        cache1.coeffs.=0.0
        exp(cache2,cache1)
    end



# power_by_squaring; slightly modified from base/intfuncs.jl
# Licensed under MIT "Expat"

     function power_by_squaring(x::Taylor0, p::Integer,cache1::Taylor0,cache2::Taylor0)
        p == 0 && return one(x,cache)
        p == 1 && return x       
        p == 2 && return square(x,cache)
        t = trailing_zeros(p) + 1
        p >>= t

        while (t -= 1) > 0
            cache1.coeffs.=0.0
            x=square(x,cache1)
        end

        cache2 = x
        while p > 0
            t = trailing_zeros(p) + 1
            p >>= t
            while (t -= 1) ≥ 0
                cache1.coeffs.=0.0
                x = square(x,cache1)
            end
            cache2 *= x
        end

        return cache2
    end
 =#

## Real power ##
function powerT(a::Taylor0, r::S,cache1::Taylor0) where {S<:Real}
    if iszero(r) #later test r=0.0
        cache1[0]=1.0
         return cache1
    end
    
  #  @show(aa)
    r == 1 && return a
    r == 2 && return square(a,cache1)
    r == 1/2 && return sqrt(a,cache1)

    l0 = findfirst(a)
   # @show l0
    lnull = trunc(Int, r*l0 )
   # @show lnull
    if (a.order-lnull < 0) || (lnull > a.order)
       # @show a.order
        return cache1  #empty
    end
    c_order = l0 == 0 ? a.order : min(a.order, trunc(Int,r*a.order))
    # @show(c_order)
    #c = Taylor0(zero(aux), c_order)
    for k = 0:c_order
        pow!(cache1, a, r, k)
    end
    # println()
    return cache1
end



# Homogeneous coefficients for real power
#= @doc doc"""
    pow!(c, a, r::Real, k::Int)

Update the `k`-th expansion coefficient `c[k]` of `c = a^r`, for
both `c` and `a` either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
c_k = \frac{1}{k a_0} \sum_{j=0}^{k-1} \big(r(k-j) -j\big)a_{k-j} c_j.
```

For `Taylor0` polynomials, a similar formula is implemented which
exploits `k_0`, the order of the first non-zero coefficient of `a`.

"""  =#
#= 
@inline function pow!(c::Taylor0{T}, a::Taylor0{T}, r::S, k::Int) where {T<:Number,S<:Real}

    if r == 0
        return one!(c, a, k)
    elseif r == 1
        return identity!(c, a, k)
    elseif r == 2
        return sqr!(c, a, k)
    elseif r == 0.5
        return sqrt!(c, a, k)
    end

    # First non-zero coefficient
    l0 = findfirst(a)
    if l0 < 0
        c[k] = zero(a[0])
        return nothing
    end

    # The first non-zero coefficient of the result; must be integer
   #=  !isinteger(r*l0) && throw(DomainError(a,
        """The 0th order Taylor0 coefficient must be non-zero
        to raise the Taylor0 polynomial to a non-integer exponent.""")) =#
    lnull = trunc(Int, r*l0 )
    kprime = k-lnull
    if (kprime < 0) || (lnull > a.order)
        @inbounds c[k] = zero(a[0])
        return nothing
    end

    # Relevant for positive integer r, to avoid round-off errors
    if isinteger(r) && r > 0 && (k > r*findlast(a))
        @inbounds c[k] = zero(a[0])
        return nothing
    end

    if k == lnull
        @inbounds c[k] = (a[l0])^r
        return nothing
    end

    # The recursion formula
    if l0+kprime ≤ a.order
        @inbounds c[k] = r * kprime * c[lnull] * a[l0+kprime]
    else
        @inbounds c[k] = zero(a[0])
    end
    for i = 1:k-lnull-1
        ((i+lnull) > a.order || (l0+kprime-i > a.order)) && continue
        aux = r*(kprime-i) - i
        @inbounds c[k] += aux * c[i+lnull] * a[l0+kprime-i]
    end
    @inbounds c[k] = c[k] / (kprime * a[l0])
    return nothing
end
 =#



#= 
## Square ##
"""
    square(a::AbstractSeries) --> typeof(a)

Return `a^2`; see [`TaylorSeries.sqr!`](@ref).
""" 
 =#

   function square(a::Taylor0,cache1::Taylor0)
        #c = Taylor0( constant_term(a)^2, a.order)
        cache1[0]=constant_term(a)^2
        for k in 1:a.order
            sqr!(cache1, a, k)
        end
        return cache1
    end




#= # Homogeneous coefficients for square
@doc doc"""
    sqr!(c, a, k::Int) --> nothing

Update the `k-th` expansion coefficient `c[k]` of `c = a^2`, for
both `c` and `a` either `Taylor0` or `TaylorN`.

The coefficients are given by

```math
\begin{eqnarray*}
c_k & = & 2 \sum_{j=0}^{(k-1)/2} a_{k-j} a_j,
    \text{ if k is odd,} \\
c_k & = & 2 \sum_{j=0}^{(k-2)/2} a_{k-j} a_j + (a_{k/2})^2,
    \text{ if k is even. }
\end{eqnarray*}
```

""" sqr! =#


    #= 
        @inline function sqr!(c::Taylor0{T}, a::Taylor0{T}, k::Int) where {T<:Number}
            if k == 0
                @inbounds c[0] = constant_term(a)^2
                return nothing
            end

            kodd = k%2
            kend = div(k - 2 + kodd, 2)
            @inbounds for i = 0:kend
                if Taylor0 == Taylor0
                    c[k] += a[i] * a[k-i]
                else
                    mul!(c[k], a[i], a[k-i])
                end
            end
            @inbounds c[k] = 2 * c[k]
            kodd == 1 && return nothing

            if Taylor0 == Taylor0
                @inbounds c[k] += a[div(k,2)]^2
            else
                sqr!(c[k], a[div(k,2)])
            end

            return nothing
        end
     =#






## Square root ##
function sqrt(a::Taylor0,cache1::Taylor0)
    l0nz = findfirst(a)
    if l0nz < 0
        cache1
    end
    cache1[0]=sqrt(a[0])
    for k = 1:a.order
        sqrt!(cache1, a, k, 0)
    end
    return cache1
end



#= @inline function sqrt!(c::Taylor0{T}, a::Taylor0{T}, k::Int, k0::Int=0) where {T<:Number}
    
    #println("k= ",k)
   # println("k0= ",k0)
    if k == k0
        @inbounds c[k] = sqrt(a[2*k0])
       println("i think this is a dead branch")
        return nothing
    end

    kodd = (k - k0)%2
   # println("kodd= ",kodd)
    kend = div(k - k0 - 2 + kodd, 2)
   # println("kend= ",kend)
    imax = min(k0+kend, a.order)
   # println("imax= ",imax)
    imin = max(k0+1, k+k0-a.order)
   # println("imin= ",imin)
    #imin ≤ imax && ( @inbounds  println("c in loop= ",c) c[k] = c[imin] * c[k+k0-imin] )
    if imin ≤ imax 
       # println("c in imin imax= ",c) 
         c[k] = c[imin] * c[k+k0-imin] 
    end
   # println("c= ",c)
    @inbounds for i = imin+1:imax
        c[k] += c[i] * c[k+k0-i]
     #   println("c in loop= ",c)
    end
    if k+k0 ≤ a.order
        @inbounds aux = a[k+k0] - 2*c[k]
      #  println("aux if <order= ",aux)
    else
        @inbounds aux = - 2*c[k]
      #  println("aux else= ",aux)
    end
    if kodd == 0
        @inbounds aux = aux - (c[kend+k0+1])^2
      #  println("aux if kodd==0  = ",aux)
    end
    @inbounds c[k] = aux / (2*c[k0])
    #println("final c inside sqrt!= ",c)
    return nothing
end
 =#

