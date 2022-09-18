
function LiQSS_integrate(::Val{O}, s::QSS_data, odep::NLODEProblem,f::Function) where {O}
  # reset_timer!()
#*********************************settings*****************************************
printCounter=[0,0]#vector{Int}
ft = s.finalTime
initTime = s.initialTime
relQ = s.dQrel
absQ = s.dQmin
#solver = s.solver
savetimeincrement=s.savetimeincrement
#*********************************qss method data*****************************************
quantum = s.quantum
#order=s.order
nextStateTime = s.nextStateTime
nextEventTime = s.nextEventTime
nextInputTime = s.nextInputTime
tx = s.tx
tq = s.tq
x = s.x
#derx=s.derx
q = s.q

t=s.t
integratorCache=s.integratorCache
taylorOpsCache=s.taylorOpsCache
cacheSize=odep.cacheSize
#=    x=Taylor0{Float64}[]
q=Taylor0{Float64}[] =#
#println("no more resize")
#*********************************problem info*****************************************
initConditions = odep.initConditions
d = odep.discreteVars
#----------to compute derivatives
jacobian = odep.jacobian
discJac = odep.discreteJacobian
#----------to compute ZC expressions
zc_jac = odep.ZC_jacobian
ZC_discJac = odep.ZC_jacDiscrete
#-----------to execute event Dependencys
evDep = odep.eventDependencies

#********************************helper values*******************************
numDiscr = length(d)
states = length(initConditions)
numberZC=size(zc_jac, 1)
numEvents=numberZC*2   #each zero crossing has 1 posEv and 1 negEv
savetime = savetimeincrement
oldsignValue = MMatrix{numberZC,2}(zeros(numberZC*2))  #usedto track if zc changed sign; each zc has a value and a sign 

#*******************************create dependencies**************************$$
SD = createDependencyMatrix(jacobian)
#SD=odep.SD
dD =createDependencyMatrix(discJac) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
SZ =createDependencyMatrix(zc_jac) 
dZ =createDependencyMatrix(ZC_discJac) # temp dependency to be used to determine HD2 and HZ2
HZ1HD1=createDependencyToEventsDiscr(dD,dZ,evDep) 
HZ2HD2=createDependencyToEventsCont(SD,SZ,evDep) 
HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
HD=unionDependency(HZ1HD1[2],HZ2HD2[2])

#= HD=odep.HD
SZ=odep.SZ
HZ=odep.HZ =#


#=   println("SD= ",SD)
println("SZ= ",SZ)
println("HZ= ",HZ)
println("HD= ",HD) =#


#@show odep.eqs
#f=@RuntimeGeneratedFunction(odep.eqs)


#@show f
zcf = Vector{Function}()
for i = 1:length(odep.zceqs)# later change to numberZC
  push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
end
#@show zcf
eventf = Vector{Function}()
for i = 1:length(odep.eventEqus)# later change to numEvents
  push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i].args[2])) 
end
#@show eventf
#*********************************  initialize          *****************************************
#= t = Taylor0(order)
integratorCache=Taylor0(zeros(order),order) #for integratestate only =#



savedVars=s.savedVars
savedTimes=s.savedTimes
  #*************************************initialize************************************
  t[0]=initTime
#  @show t
#there is a taylor creation after this
#for k = 1:order # compute initial derivatives for x and q (similar to a recursive way )
  #for i = 1:states
     # q computed from x and it is going to be used in the next x
  #end

  for i = 1:states
    #q[i].coeffs .= x[i].coeffs
    q[i][0] = x[i][0]
   # q[i][1] = x[i][1]
   # q[i][2] = x[i][2]
  end
  for i = 1:states
    clearCache(taylorOpsCache,cacheSize)
    f(i,q,d, t ,taylorOpsCache)
    x[i][1]=taylorOpsCache[1][0]

   # x[i][1]=f(i,q,d, t ,taylorOpsCache)[0]
  end
  for i = 1:states
    q[i][1] = x[i][1]
  end
  for i = 1:states
    clearCache(taylorOpsCache,cacheSize)
    f(i,q,d, t ,taylorOpsCache)
    x[i][2]=taylorOpsCache[1][1]/2
    
    #x[i][2]=f(i,q,d, t ,taylorOpsCache)[1]/2
   # q[i][2] = x[i][2]# not needed q is kept the same order to avoid resize but derderQ is never updated and it  is always 0
   # x[i].coeffs[k+1] = ((differentiate( derx[i], k - 1)).coeffs[1]) / factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
   #=  f[i](integratorCache,q,d, t + initTime,taylorOpsCache)
    x[i].coeffs[k+1] = ((differentiate(integratorCache, k - 1)).coeffs[1]) / factorial(k) =#
  end

for i = 1:states

  #println("before bcast who is resizing?")
    savedVars[i][1].coeffs .= x[i].coeffs  #to be changed  1 2 3 ?
    #= savedVars[i][1].coeffs[1] = x[i].coeffs[1]
    savedVars[i][1].coeffs[2] = x[i].coeffs[2]
    savedVars[i][1].coeffs[3] = x[i].coeffs[3] =#

  quantum[i] = relQ * abs(x[i].coeffs[1]) 
  if quantum[i] < absQ
    quantum[i] = absQ
  end
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  clearCache(taylorOpsCache,cacheSize)

  f(i,q,d,t,taylorOpsCache) #+t alloc   change to addT
  computeNextInputTime(Val(O), i, initTime, 0.0,taylorOpsCache[1], nextInputTime, x,  quantum)
  #+t alloc   change to addT
  computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)

end
#println("before zc who is resizing?")
for i=1:numberZC
  clearCache(taylorOpsCache,cacheSize)
  output=zcf[i](x,d,t,taylorOpsCache).coeffs[1] #test this evaluation
  oldsignValue[i,2]=output #value
  oldsignValue[i,1]=sign(output) #sign modify 
  computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum,printCounter)
end

#= 
println("before while situation")
@show x
@show q
@show quantum
@show nextStateTime 
@show nextEventTime 
@show nextInputTime =#


###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime
count = 1 # not zero because intial value took 0th position
len=length(savedTimes)
#@show len
#printcount=0
while simt < ft #&& printcount < 3
  #printcount += 1
  sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2]
  index = sch[1]
  t[0]=simt
 # println("index= ",index)
  ##########################################state########################################
  if sch[3] == :ST_STATE
  
  #=   if printCounter[2] > 0 && 1.4386 < simt < 1.4388
      println("at time ",simt," state changed after the event is x", index)

      println("nexttime= ", nextStateTime)
      println("nexteventtime= ", nextEventTime)
      println("nextinputtime= ", nextInputTime)
      println("x= ",x)
      
    end =#
    elapsed = simt - tx[index]
   integrateState(Val(O),x[index],integratorCache,elapsed)
    quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
    if quantum[index] < absQ
      quantum[index] = absQ
    end
 
    tx[index] = simt
    for k = 1:O
      q[index].coeffs[k] = x[index].coeffs[k] #updateQ
    end
    tq[index] = simt #tq needed for higher orders
     computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
  #=    println("state: x= ",x)
     println("state: q= ",q)
     println("state: nexttime= ",nextStateTime)
     println("quantum= ",quantum) =#
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0           
        elapsed = simt - tx[j]
        if elapsed > 0
          #"evaluate" x at new time only...derivatives get updated next using computeDerivative()
          x[j].coeffs[1] = x[j](elapsed)
          q[j].coeffs[1] = x[j].coeffs[1]
          tx[j] = simt
          tq[j] = simt
        end
       # @timeit "state-clearCache" clearCache(taylorOpsCache,cacheSize)
         #normally and later i should update other q (integrate q=q+e derQ  for higher orders)
        #integrate!(x[j],f(j,q,d,t,taylorOpsCache),x[j].coeffs[1])
       # @timeit "state-f"  
        #f(j,q,d,t,taylorOpsCache)
        #@timeit "state-computeder" 
        f(j,q,d,t,taylorOpsCache)
        computeDerivative2(Val(O), x[j], taylorOpsCache[1])
         reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
        #=  if printCounter[2] > 0
          printCounter[2]-=1
         println("x[$j]= ",x[j])
         end =#
      #=    println("state: after compder x= ",x)
         println("state: q= ",q)
         println("state: nexttime= ",nextStateTime) =#
      end#end if j!=0
    end#end for SD
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache,cacheSize)
      #  println("zcf[j](x,d,t,taylorOpsCache)= ",zcf[j](x,d,t,taylorOpsCache))
       # computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum,printCounter)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum,printCounter)
     #   println("state: nexteventtime= ",nextEventTime)
      end  #end if j!=0
    end#end for SZ
    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself  
   # println("begin input:who is resizing?")
    elapsed = simt - tx[index]    
    clearCache(taylorOpsCache,cacheSize)
    
    f(index,q,d,t,taylorOpsCache)
     computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
     
    # computeNextInputTime(Val(2), index, simt, elapsed,f(index,q,d,t,taylorOpsCache) , nextInputTime, x,  quantum)
   #  println("middle input:who is resizing?")
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0             
        elapsed = simt - tx[j]
        if elapsed > 0

          x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
          tx[j] = simt
        end
        quantum[j] = relQ * abs(x[j].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
        if quantum[j] < absQ
          quantum[j] = absQ
        end
        clearCache(taylorOpsCache,cacheSize)
        
         #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        # f(j,q,d,t,taylorOpsCache)
       # integrate!(x[j][j],x[j].coeffs[1]) 
       f(j,q,d,t,taylorOpsCache)
        computeDerivative2(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end if j!=0
    end#end for
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache,cacheSize)
        
        computeNextEventTime(j,zcf[j](q,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum,printCounter) 
      end  
     # println("end input:who is resizing?")
    end
  #################################################################event########################################
  else
    #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
    modifiedIndex=0
   # clearCache(taylorOpsCache,cacheSize)#     *****************
 #=   if printCounter[1]>0 && 1.4386 < simt < 1.4388
    println("event occured at time= ",simt)
    println("zcf=",zcf[index](x,d,t,taylorOpsCache))
    println("q1=   ",q[1])
    println("d= ",d)
    println("index= ",index)
   end =#
    if (zcf[index](x,d,t,taylorOpsCache).coeffs[1])>0       # sign is not needed here
      modifiedIndex=2*index-1   # the  event that just occured is at  this index
    else
      modifiedIndex=2*index
    end    
       
    #=
      since i am checking if zc>0 then posEV Funct1 else negEvFunct2....it is better to store the whole thing as an event function
      that way i can have more complicated events without worrying about parsing...one thing: is it cheaper to store one giant event in a function
      or too many small events...i think one giant event is cheaper.
    =#
   # clearCache(taylorOpsCache,cacheSize)
    
    eventf[modifiedIndex](q,d,t,taylorOpsCache)
#=     if printCounter[1]>0 
     
      println("modified index= ",modifiedIndex)
      println("d= ",d)
      
      printCounter[2] =1
  
     end =# 
   # println("modified index= ",modifiedIndex)
    #println("eventf ",eventf)
    #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
    nextEventTime[index]=Inf   #not tested
   # printCounter[1]=3
    for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
        j = HD[modifiedIndex][i] # 
        if j != 0
         #=  if printCounter[1]>0 && 1.4386 < simt < 1.4388
            #println("event occured")
            println("the state var affected by event is x$j")
        
        
           end =#   
                elapsed = simt - tx[j]
              
                if elapsed > 0  # if event triggere by change of sign and time=now then elapsed=0

                   #=  if printCounter[1]>0 
                      #println("event occured")
                      println("elapsed appear?= ",elapsed) #only appear when zc depends on x1 and it affects derx2 && derx1 doesn't depend on x1
                  
                  
                     end  =#
                   
                  x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                  q[j].coeffs[1] = x[j].coeffs[1]
                  tx[j] = simt
                end
                #computeDerivative(Val(2), j,x, f(j,q(simt + t),d,simt + t)) #updating q as a vector updates all components q[0] q[1]...
                #but some q is already updated!!!causes error. if q[1] updated at 0.5 then q(t+0.5) makes another change
               #=  quantum[j] = relQ * abs(x[j].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
                if quantum[j] < absQ
                  quantum[j] = absQ
                end =#
                clearCache(taylorOpsCache,cacheSize)
             # f(j,q,d,t,taylorOpsCache)
            #  integrate!(x[j][j],x[j].coeffs[1])
           #=  if printCounter[1]>0 && 1.4386 < simt < 1.4388
              println("x$j before der= ",x)   
              println("q$j before der= ",q)
              println(" f(j,q,d,t,taylorOpsCache)= ", f(j,q,d,t,taylorOpsCache))
              clearCache(taylorOpsCache,cacheSize)
              println("quantum= ",quantum)
              println("nextStateTime= ",nextStateTime)
              println("simt= ",simt)
              end =#
              f(j,q,d,t,taylorOpsCache)
              computeDerivative2(Val(O), x[j], taylorOpsCache[1])
              reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)  
             #=  if printCounter[1]>0  && 1.4386 < simt < 1.4388
              println("x$j after event= ",x)   
              printCounter[1]-=1 
              println("nexttime= ",nextStateTime)
              end =#
        end     
    end
    for i = 1:length(HZ[modifiedIndex])
          j = HZ[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
            if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)
            
            computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum,printCounter)
          end        
    end
   
  end#end state input event
  if simt > savetime
      count += 1
      savetime += savetimeincrement #next savetime
      if len<count
        len=count*2
        for i=1:states
          resize!(savedVars[i],len)
        end
        resize!(savedTimes,len)
      end
      for k = 1:states
         #=  temp = Taylor0(interpolationOrder)
          for j = 1:interpolationOrder
            temp.coeffs[j] = x[k].coeffs[j]
          end
          # println("temp= ",temp)
          push!(savedVars[k], temp) =#
          #push!(savedVars[k], x[k][0])
         #=  if length(savedVars[k])<count
            resize!(savedVars[k],count*2)
          end =#
          #savedVars[k][count]=x[k][0]
          savedVars[k][count].coeffs .=x[k].coeffs
         #=  savedVars[k][count][0]=x[k][0]
          savedVars[k][count][1]=x[k][1]
          savedVars[k][count][2]=x[k][2] =#
    
      end
     # push!(savedTimes, simt)
     savedTimes[count]=simt
  end#end if save
end#end while
#print_timer()
#= @show len
@show count =#
for i=1:states# throw away empty points
  resize!(savedVars[i],count)
end
resize!(savedTimes,count)
(savedTimes, savedVars)
end#end integrate




#= function f(j::Int,q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
  if j == 1
   createT(q[2],cache[1])
  return nothing
else
    (subT(negateT(q[1], cache[2]), q[2], cache[1]))#.coeffs
   return nothing
end
end =#

#struct Solution
