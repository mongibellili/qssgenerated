
 function updateScheduler(nextStateTime::MVector{T,Float64},nextEventTime :: MVector{Z,Float64},nextInputTime :: Vector{Float64} )where{T,Z}   #later MVect for nextInput
    minStateTime=Inf
    minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
    minEventTime=Inf
    minEvent_index=0
    minInputTime=Inf
    minInput_index=0

    returnedVar=() #  used to print something if something is bad
    for i=1:T
        if nextStateTime[i]<minStateTime
            minStateTime=nextStateTime[i]
            minState_index=i
        end
    end
    for i=1:Z
        if nextEventTime[i] < minEventTime
            minEventTime=nextEventTime[i]
            minEvent_index=i
        end
    end
    for i=1:T
        if nextInputTime[i] < minInputTime
            minInputTime=nextInputTime[i]
            minInput_index=i
        end
    end

    if minEventTime<minStateTime
       if minInputTime<minEventTime
         returnedVar=(minInput_index,minInputTime,:ST_INPUT)
       else
          returnedVar=(minEvent_index,minEventTime,:ST_EVENT)     
       end
    else      
        if minInputTime<minStateTime
            returnedVar=(minInput_index,minInputTime,:ST_INPUT)
        else
             returnedVar=(minState_index,minStateTime,:ST_STATE)
        end
    end
    if returnedVar[1]==0
        println("scheduler *******nextEventTime= ",nextEventTime)
        println("scheduler *******nextInputTime= ",nextInputTime)
        println("scheduler *******nextstateTime= ",nextStateTime)
        println("also the whole system may be static! developer: fill remaing points with same values and exit!")
    end
    return returnedVar 
end
