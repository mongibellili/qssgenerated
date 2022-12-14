@inline function smallest(args::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	w = v
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	d = zero(F)
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	if low < zero(F)
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, cohigh, coeff), muladd(high, colow, coeff)
			elseif cohigh < zero(F)
				muladd(high, cohigh, coeff), muladd(low, colow, coeff)
			else
				muladd(low, cohigh, coeff), muladd(low, colow, coeff)
			end
		end
	else
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, colow, coeff), muladd(high, cohigh, coeff)
			elseif cohigh < zero(F)
				muladd(high, colow, coeff), muladd(low, cohigh, coeff)
			else
				muladd(high, colow, coeff), muladd(high, cohigh, coeff)
			end
		end	
	end
    colow, cohigh
end

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > zero(F)
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < zero(F)
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		return typemax(F)
	elseif N == 2
		@inbounds ret = -coeffs[2] / coeffs[1]
		if ret < zero(F)
			return typemax(F)
		else
			return ret
		end
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly??? = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	MM = smallest(poly)
	if MM > zero(F) return typemax(F) end
	domlow, domhigh = zero(F), one(F) - MM
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom???low, codom???high = posintervalhorner(domlow, domhigh, poly???)
		#@show comid codom???low codom???high
		if codom???low < zero(F) < codom???high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom???low, mid - comid / codom???high, domhigh
			else
				domlow, mid - comid / codom???high, mid - comid / codom???low, domhigh
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
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid??? = hornerd(mid, poly)
					delta = comid / comid???
					newmid = mid - delta
					if abs(delta) < 1.0e-8mid 
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
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
	return typemax(F)
end

function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::Ptr{F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		return 0
	elseif N == 2
		@inbounds res[1] = -coeffs[2] / coeffs[1]
		return 1
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly??? = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	mm = zero(F)
	MM = zero(F)
	s = -one(F)
	for i in 2:N
		@inbounds coeff = poly[i]
		if coeff < MM; MM = coeff end
		_coeff = s * coeff
		if _coeff < mm; mm = _coeff end
		s = -s
	end
	if mm == zero(F) && MM == zero(F) return 0 end
	index = 0
	domlow, domhigh = if mm < zero(F)
		if MM < zero(F)
			index = 1
			unsafe_store!(doms, (zero(F), one(F) - MM), 1)
		end
		mm - one(F), zero(F)
	else
		zero(F), one(F) - MM
	end
    counter = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom???low, codom???high = intervalhorner(domlow, domhigh, poly???)
		#@show mid comid codom???low codom???high
		if codom???low < zero(F) < codom???high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom???low, mid - comid / codom???high, domhigh
			else
				domlow, mid - comid / codom???high, mid - comid / codom???low, domhigh
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
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid??? = hornerd(mid, poly)
					delta = comid / comid???
                    #@show mid comid comid???
					newmid = mid - delta
					if abs(delta) < 1.0e-8abs(mid) 
						counter += 1
                        unsafe_store!(res, newmid, counter)
                        break
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 || counter == N-1 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return counter
end
#coeffs=NTuple{4,Float64}((1.33,12.65,1.0,-6.0))

pp=pointer(Vector{NTuple{2,Float64}}(undef, 3))
#@show pp
#@btime gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#64.912 ns (1 allocation: 32 bytes)
#@show gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#2.414213562506035#@show quadraticsmallestpositiveroot(1.0, -2.0, -1.0)#1.649 ns (0 allocations: 0 bytes)
@show smallestpositiverootintervalnewtonregulafalsi((1.0,-4.0,4.0),pp)