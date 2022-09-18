# This file is part of the Taylor0Series.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#
#= 
"""
    AbstractSeries{T<:Number} <: Number

Parameterized abstract type for [`Taylor0`](@ref),
[`HomogeneousPolynomial`](@ref) and [`TaylorN`](@ref).
""" =#
abstract type AbstractSeries{T<:Number} <: Number end


## Constructors ##

######################### Taylor0
#= """
    Taylor0{T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in one independent variable.

**Fields:**

- `coeffs :: Array{T,1}` Expansion coefficients; the ``i``-th
    component is the coefficient of degree ``i-1`` of the expansion.
- `order  :: Int` Maximum order (degree) of the polynomial.

Note that `Taylor0` variables are callable. For more information, see
[`evaluate`](@ref).
""" =#
struct Taylor0{T<:Number} <: AbstractSeries{T}
    coeffs :: Array{T,1}
    order :: Int

    ## Inner constructor ##
    function Taylor0{T}(coeffs::Array{T,1}, order::Int) where T<:Number
        #resize_coeffs1!(coeffs, order)
        return new{T}(coeffs, order)
    end
end

## Outer constructors ##
Taylor0(x::Taylor0{T}) where {T<:Number} = x
Taylor0(coeffs::Array{T,1}, order::Int) where {T<:Number} = Taylor0{T}(coeffs, order)
Taylor0(coeffs::Array{T,1}) where {T<:Number} = Taylor0(coeffs, length(coeffs)-1)
function Taylor0(x::T, order::Int) where {T<:Number}
    v = fill(zero(x), order+1)
    v[1] = x
    return Taylor0(v, order)
end

# Methods using 1-d views to create Taylor0's
Taylor0(a::SubArray{T,1}, order::Int) where {T<:Number} = Taylor0(a.parent[a.indices...], order)
Taylor0(a::SubArray{T,1}) where {T<:Number} = Taylor0(a.parent[a.indices...])


# Shortcut to define Taylor0 independent variables
#= """
    Taylor0([T::Type=Float64], order::Int)

Shortcut to define the independent variable of a `Taylor0{T}` polynomial of
given `order`. The default type for `T` is `Float64`.

```julia
julia> Taylor0(16)
 1.0 t + 𝒪(t¹⁷)

julia> Taylor0(Rational{Int}, 4)
 1//1 t + 𝒪(t⁵)
```
""" =#
Taylor0(::Type{T}, order::Int) where {T<:Number} = Taylor0( [zero(T), one(T)], order)
Taylor0(order::Int) = Taylor0(Float64, order)



# A `Number` which is not an `AbstractSeries`
const NumberNotSeries = Union{Real,Complex}

# A `Number` which is not `TaylorN` nor a `HomogeneousPolynomial`
const NumberNotSeriesN = Union{Real,Complex,Taylor0}

## Additional Taylor0 outer constructor ##
Taylor0{T}(x::S) where {T<:Number,S<:NumberNotSeries} = Taylor0([convert(T,x)], 0)
