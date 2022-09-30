using BenchmarkTools


@inline function smallest(args::NTuple{4,Float64})::Float64# where {N, F <: AbstractFloat}
	@inbounds m = args[1]
	for i in 2:4
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end
#= @generated function horner(x::Float64, coeff1::Float64, coeffs::Float64...)::Float64# where F <: AbstractFloat
    l = length(coeffs)
   l === 0 && return quote coeff1 end
   ex = :(coeff1)
   for i in 1:l
       ex = :(@inbounds muladd(x, $ex, coeffs[$i]))
   end
   ex |> flatten
end
@inline function horner(x::Float64, coeff1::Float64,coeff2::Float64,coeff3::Float64,coeff4::Float64 )# where {N, F <: AbstractFloat}
    v = muladd(x, muladd(x, muladd(x, coeff1, coeff2), coeff3), coeff4)
end =#

#= @generated function horner(x::Float64, coeffs::NTuple{4,Float64})::Float64# where {N, F <: AbstractFloat}
	ex = :(coeffs[1])
	for i in 2:4
		ex = :(@inbounds muladd(x, $ex, coeffs[$i]))
	end
	ex 
end =#
@inline function horner(x::Float64, coeffs::NTuple{4,Float64})::Float64# where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	for i in 2:4
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::Float64, y::Float64, coeffs::NTuple{4,Float64})::Tuple{Float64, Float64}# where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	w = v
	for i in 2:4
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::Float64, coeffs::NTuple{4,Float64})::Tuple{Float64, Float64}# where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	d = 0.0
	for i in 2:4
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function posintervalhorner(low::Float64, high::Float64, coeffs::NTuple{3,Float64})::Tuple{Float64, Float64}# where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	for i in 2:3
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > 0.0
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < 0.0
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end


#= function lagrangeinversequadratic(x1::Float64,x2::Float64,x3::Float64,p1::Float64, p2::Float64, p3::Float64)::Float64
    p12=p1-p2
    p13=p1-p3
    p23=p2-p3
    finv=(p2*p3*x1*p23-p1*p3*x2*p13+p1*p2*x3*p12)/(p12*p13*p23)
end =#

#= coeffs=NTuple{4,Float64}((-1.33,2.0,1.0,6.0))
@inbounds _coeff1 = inv(coeffs[1])
poly = ntuple(4) do i
    @inbounds _coeff1 * coeffs[i]
end

poly′ = ntuple(3) do i
    @inbounds (3) * poly[i]
end =#
#@btime smallest(poly)#3.307 ns (0 allocations: 0 bytes)
#= domlow=0.0
domhigh=0.75
mid=0.37

coeff1=1.0
coeff2=2.0/-1.33
coeff3=1.0/-1.33
coeff4=6.0/-1.33 =#
#@btime horner(mid, 1.0,-1.50375939849,-0.7518796992,-4.511278195)#5.040 ns (0 allocations: 0 bytes)
#@btime horner(mid, coeff1,coeff2,coeff3,coeff4)#5.040 ns (0 allocations: 0 bytes)
#@btime posintervalhorner(domlow, domhigh, poly′)#6.527 ns (0 allocations: 0 bytes)
#@btime horner2(domlow, domhigh, poly)#5.717 ns (0 allocations: 0 bytes)
#@btime hornerd(mid, poly)#5.462 ns (0 allocations: 0 bytes)
#@btime horner(mid,(1.0, -1.5037593984, -0.751879699, -4.511278195))#5.040 ns (0 allocations: 0 bytes)
function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{4,Float64}, doms::Ptr{NTuple{2,Float64}})::Float64# where {N, F <: AbstractFloat}
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(4) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(3) do i
        @inbounds (3) * poly[i]
    end
	MM = smallest(poly)
	if MM > 0.0 return Inf end
	domlow, domhigh = 0.0, 1.0- MM
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
		#@show comid codom′low codom′high
		if codom′low < 0.0 < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid < 0.0
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					unsafe_store!(doms, (rightlow, righthigh), index)
				end
				domlow, domhigh = leftlow, lefthigh
				continue
			elseif rightlow < righthigh
				domlow, domhigh = rightlow, righthigh
				continue
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
			#@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < 0.0
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
					newmid = mid - delta
					newmidy=horner(newmid,poly)
					if abs(delta) < 1.0e-6mid 
						#newmid=lagrangeinversequadratic(domlow, newmid,mid,codomlow, newmidy,comid)
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < 0.0
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return Inf
end


#coeffs=NTuple{4,Float64}((-1.33,2.0,1.0,6.0))
coeffs=NTuple{4,Float64}((-1.33,2.0,1.0,6.0))

pp=pointer(Vector{NTuple{2,Float64}}(undef, 3))
#pp=pointer(Vector{NTuple{2,Float64}}(undef, 1))#signal (11): Segmentation fault
#@show pp
#@btime gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#64.912 ns (1 allocation: 32 bytes)
#@show gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#2.414213562506035#@show quadraticsmallestpositiveroot(1.0, -2.0, -1.0)#1.649 ns (0 allocations: 0 bytes)
#@show smallestpositiverootintervalnewtonregulafalsi((-1.33,2.0,1.0,6.0),pp)
#@btime smallestpositiverootintervalnewtonregulafalsi((-1.33,2.0,1.0,6.0),pp)# 93.687 ns (0 allocations: 0 bytes)
@btime smallestpositiverootintervalnewtonregulafalsi(coeffs,pp)# 93.687 ns (0 allocations: 0 bytes)
#@btime smallestpositiverootintervalnewtonregulafalsi((-1.33,2.0,1.0,6.0),pp)#126.486 ns (2 allocations: 64 bytes)....2.5155562382541037...cubic:18ns
#@btime gensmallestpositiverootintervalnewton(1.33,2.0,-2.63,6.0)#32.639 ns (0 allocations: 0 bytes) ....Inf
#@show gensmallestpositiverootintervalnewton(2.0,-3.0,-0.5,-5.0)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0,2.2)#171.461 ns (2 allocations: 64 bytes)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0,2.2,3.14)#220.877 ns (2 allocations: 64 bytes)
