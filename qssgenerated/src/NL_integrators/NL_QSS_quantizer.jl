@inline function integrateState(::Val{1}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
  x.coeffs[1] = x(elapsed)
end
@inline function integrateState(::Val{2}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
  x.coeffs[1] = x(elapsed)
  differentiate!(cacheT,x)
  x.coeffs[2] = cacheT(elapsed)
  #cacheT.coeffs.=0.0 #clear cache
end


function computeDerivative( ::Val{1}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64  )
   x.coeffs[2] =f(elap)
   #x.coeffs[2] =f.coeffs[1]
   return nothing
end
#= function computeDerivative( ::Val{2}  ,x::Taylor0{Float64},f::Taylor0{Float64}  ) # approx deriv ...same alloc gain in 11us for FT=5.0  58---69
   x.coeffs[2] =f.coeffs[1]
   x.coeffs[3] =f.coeffs[2]/2
   return nothing
end =#
function computeDerivative( ::Val{2}  ,x::Taylor0{Float64},f::Taylor0{Float64},cache::Taylor0{Float64},elap::Float64 )
  x.coeffs[2] =f(elap)
  differentiate!(cache,f)
  x.coeffs[3]=cache(elap)/2
  return nothing
end


function computeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
    if (x[i].coeffs[2]) != 0
        tempTime=max((abs(quantum[i] / ((x[i].coeffs[2])))),absDeltaT)
        if tempTime!=absDeltaT #normal
            nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
        else#usual (quant/der) is very small
          x[i].coeffs[2]=sign(x[i].coeffs[2])*(abs(quantum[i])/absDeltaT)# adjust  derivative if it is too high
          nextTime[i] = currentTime + tempTime
          println("smalldelta")
        end
    else
      nextTime[i] = Inf
    end
    return nothing
end

function computeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    absDeltaT=1e-12 # minimum deltaT to protect against der=Inf coming from sqrt(0) for example...similar to min ΔQ
      if (x[i].coeffs[3]) != 0
          tempTime=max(sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))),absDeltaT)
          if tempTime!=absDeltaT #normal
              nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
          else#usual sqrt(quant/der) is very small
            x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
            nextTime[i] = currentTime + tempTime
          end
      else
        nextTime[i] = Inf
      end
      return nothing
end


function reComputeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], -x[index].coeffs[2]]
  time1 = currentTime + minPosRoot(coef, Val(1))
  coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  time2 = currentTime + minPosRoot(coef, Val(1))
  nextTime[index] = time1 < time2 ? time1 : time2
  #later guard against very small Δt like in computenext
  return nothing
end

function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
  coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])*2]
  time1 = currentTime + minPosRoot(coef, Val(2))
  coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
  time2 = currentTime + minPosRoot(coef, Val(2))
  nextTime[index] = time1 < time2 ? time1 : time2
  #later guard against very small Δt like in computenext
  return nothing
end
  
function computeNextInputTime(::Val{1}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T<:Int64}
    df=0.0
    oldDerX=x[i].coeffs[2]
    newDerX=tt.coeffs[1] 
    if elapsed > 0
      df=(newDerX-oldDerX)/elapsed
    else
      df= quantum[i]*1e6
    end     
    if df!=0.0
       nextInputTime[i]=currentTime+sqrt(abs(quantum[i] / df))
    else
      nextInputTime[i] = Inf
      println("nextinput inf")
    end
    return nothing
end


  
function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T<:Int64}
  df=0.0
      oldDerDerX=((x[i].coeffs[3])*2.0)
    newDerDerX=(tt.coeffs[2])
      if elapsed > 0.0
        df=(newDerDerX-oldDerDerX)/(elapsed/2)
      else
        df= quantum[i]*1e12
      end       
     if df!=0.0
      nextInputTime[i]=currentTime+cbrt((abs(quantum[i]/df)))     
     else
      nextInputTime[i] = Inf
      end
      return nothing
end
function computeNextEventTime(j::Int,ZCFun::Taylor0{Float64},oldsignValue,currentTime,  nextEventTime, quantum::Vector{Float64},printCounter::Vector{Int}) #later specify args
  if oldsignValue[j,1] != sign(ZCFun.coeffs[1])
    nextEventTime[j]=currentTime 
  else
    #nextEventTime[j] =currentTime + minPosRoot(ZCFun.coeffs, Val(2)) #Inf  # we can estimate the time. this is important when zc depends only on time   
    nextEventTime[j]=Inf
  end
  oldsignValue[j,1]=sign(ZCFun.coeffs[1])#update the values
  oldsignValue[j,2]=ZCFun.coeffs[1]
end

