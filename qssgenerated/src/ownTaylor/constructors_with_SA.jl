# This file is part of the TaylorSASeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
using StaticArrays
abstract type AbstractSeries <: Number end
## Constructors ##
######################### TaylorSA
#= 
- `coeffs :: MVector{1,Float64}` Expansion coefficients; the ``i``-th
    component is the coefficient of degree ``i-1`` of the expansion.
- `order  :: Int` Maximum order (degree) of the polynomial.
""" =#
struct TaylorSA{V} #<: AbstractSeries
    coeffs :: MVector{V,Float64}
    order :: Int
    function TaylorSA{Float64}(coeffs::MVector{V,Float64}, order::Int) where {V}
        #resize_coeffs1!(coeffs, order)
        return new{V}(coeffs, order)
    end
end
TaylorSA(x::TaylorSA{Float64})  = x
TaylorSA(coeffs::MVector{V,Float64}, order::Int) where {V} = TaylorSA{Float64}(coeffs, order)
TaylorSA(coeffs::MVector{V,Float64}) where {V} = TaylorSA(coeffs, length(coeffs)-1)
function TaylorSA(x::Float64, order::Int) 
    #MVector{order+1,Int}(vect)
    v=@MVector zeros(order+1)
    #v = fill(zero(x), order+1)
    v[1] = x
    return TaylorSA(v, order)
end
TaylorSA(order::Int) = TaylorSA(1.0, order)
vect=(1.2,2.2,3.0)
coeffs=MVector{3,Float64}(vect)
t1= TaylorSA(coeffs,2)
t2= TaylorSA(2.3,5)
t3= TaylorSA(5)
#@show t1, t2, t3
import Base.:+
function (+)(a::TaylorSA{T}, b::T) where {T}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (+)(a[0], b)
    return TaylorSA(coeffs, a.order)
end
t4=t1+12.0
@show t4
