# using TimerOutputs
function QSS_integrate(::Val{O}, s::QSS_data{T,Z}, odep::NLODEProblem{T,D,Z,Y},f::Function) where {O,T,D,Z,Y}
#  reset_timer!()
#*********************************settings*****************************************
#printCounter=[0,0]#vector{Int} fort debugging to be deleted
ft = s.finalTime
initTime = s.initialTime
relQ = s.dQrel
absQ = s.dQmin
savetimeincrement=s.savetimeincrement
#*********************************qss method data*****************************************
quantum = s.quantum
nextStateTime = s.nextStateTime
nextEventTime = s.nextEventTime
nextInputTime = s.nextInputTime
tx = s.tx
tq = s.tq
x = s.x
q = s.q
t=s.t
savedVars=s.savedVars
savedTimes=s.savedTimes
integratorCache=s.integratorCache
taylorOpsCache=s.taylorOpsCache
cacheSize=odep.cacheSize
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
#states = length(initConditions)
#numberZC=size(zc_jac, 1)
numEvents=Z*2   #each zero crossing has 1 posEv and 1 negEv
savetime = savetimeincrement
oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 
#*******************************create dependencies**************************$$
SD = createDependencyMatrix(jacobian)
dD =createDependencyMatrix(discJac) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
SZ =createDependencyMatrix(zc_jac) 
dZ =createDependencyMatrix(ZC_discJac) # temp dependency to be used to determine HD2 and HZ2
HZ1HD1=createDependencyToEventsDiscr(dD,dZ,evDep) 
HZ2HD2=createDependencyToEventsCont(SD,SZ,evDep) 
HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
HD=unionDependency(HZ1HD1[2],HZ2HD2[2])
#=   println("SD= ",SD)
println("SZ= ",SZ)
println("HZ= ",HZ)
println("HD= ",HD) =#
zcf = Vector{Function}()
for i = 1:length(odep.zceqs)# later change to numberZC
  push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
end
eventf = Vector{Function}()
for i = 1:length(odep.eventEqus)# later change to numEvents
  push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i].args[2])) 
end

#######################################compute initial values##################################################
n=1
for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
  n=n*k
   for i = 1:T
      q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
    end
    for i = 1:T
      clearCache(taylorOpsCache,cacheSize)
       f(i,q,d, t ,taylorOpsCache)
       (ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1))
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
end

for i = 1:T
  savedVars[i][1].coeffs .= x[i].coeffs  #to be changed  1 2 3 ?
  quantum[i] = relQ * abs(x[i].coeffs[1]) 
  if quantum[i] < absQ
    quantum[i] = absQ
  end
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  clearCache(taylorOpsCache,cacheSize)
  f(i,q,d,t,taylorOpsCache) #+t alloc   change to addT
  computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)
end

for i=1:Z
  clearCache(taylorOpsCache,cacheSize)
  output=zcf[i](x,d,t,taylorOpsCache).coeffs[1] #test this evaluation
  oldsignValue[i,2]=output #value
  oldsignValue[i,1]=sign(output) #sign modify 
  computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum)#,printCounter)
end

###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime
count = 1 # not zero because intial value took 0th position
len=length(savedTimes)
#printcount=0
while simt < ft #&& printcount < 10
  #printcount+=1
  sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2]
  index = sch[1]
  t[0]=simt
  ##########################################state########################################
  if sch[3] == :ST_STATE
    elapsed = simt - tx[index]
    #@timeit "state-integrateState" 
    integrateState(Val(O),x[index],integratorCache,elapsed)
    quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
    if quantum[index] < absQ
      quantum[index] = absQ
    end
    tx[index] = simt
    #@timeit "updateQ" 
    for k = 1:O
      q[index].coeffs[k] = x[index].coeffs[k] #updateQ
    end
    tq[index] = simt #tq needed for higher orders down before computing f(q,d,t)
    #@timeit "state-computenext"
     computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0           
        elapsed = simt - tx[j]
        if elapsed > 0
          #"evaluate" x only at new time ...derivatives get updated next using computeDerivativ()
          x[j].coeffs[1] = x[j](elapsed)
          q[j].coeffs[1] = q[j](elapsed)
          tx[j] = simt
          tq[j] = simt
        end
        f(j,q,d,t,taylorOpsCache)
        #@timeit "state-compder" 
         computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
        #computeDerivative(Val(O), x[j], taylorOpsCache[1])
       # @timeit "state-recomputenext"  
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end if j!=0
    end#end for SD
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache,cacheSize)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
      end  #end if j!=0
    end#end for SZ
    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself  

    elapsed = simt - tx[index]    
    clearCache(taylorOpsCache,cacheSize)   
    f(index,q,d,t,taylorOpsCache)
     computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0             
        elapsed = simt - tx[j]
        if elapsed > 0

          x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivati()
          tx[j] = simt
        end
        quantum[j] = relQ * abs(x[j].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
        if quantum[j] < absQ
          quantum[j] = absQ
        end
        clearCache(taylorOpsCache,cacheSize)
        f(j,q,d,t,taylorOpsCache)
        computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
       # computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end if j!=0
    end#end for
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache,cacheSize)
        
        computeNextEventTime(j,zcf[j](q,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter) 
      end  
     # println("end input:who is resizing?")
    end
  #################################################################event########################################
  else
    #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 

    modifiedIndex=0
    if (zcf[index](x,d,t,taylorOpsCache).coeffs[1])>0       # sign is not needed here
      modifiedIndex=2*index-1   # the  event that just occured is at  this index
    else
      modifiedIndex=2*index
    end       
    eventf[modifiedIndex](q,d,t,taylorOpsCache)
    #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
    nextEventTime[index]=Inf   #investigate more
    for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
        j = HD[modifiedIndex][i] # 
        if j != 0
              elapsed = simt - tx[j]             
              if elapsed > 0  # if event triggere by change of sign and time=now then elapsed=0
                  # println("elapsed appear?= ",elapsed) #only appear when zc depends on x1 and it affects derx2 && derx1 doesn't depend on x1             
                  x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivativ()
                  q[j].coeffs[1] = x[j].coeffs[1]
                  tx[j] = simt
              end
                clearCache(taylorOpsCache,cacheSize)

              f(j,q,d,t,taylorOpsCache)
              computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
             # computeDerivative(Val(O), x[j], taylorOpsCache[1])
              reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)  
        end     
    end
    for i = 1:length(HZ[modifiedIndex])
          j = HZ[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
            if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)           
            computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
          end        
    end
  end#end state/input/event
  if simt > savetime
      count += 1
      savetime += savetimeincrement #next savetime
      if len<count
        len=count*2
        for i=1:T
          resize!(savedVars[i],len)
        end
        resize!(savedTimes,len)
      end
      for k = 1:T
          savedVars[k][count].coeffs .=x[k].coeffs 
      end
     savedTimes[count]=simt
  end#end if save
end#end while
for i=1:T# throw away empty points
  resize!(savedVars[i],count)
end
# print_timer()
resize!(savedTimes,count)
Sol(savedTimes, savedVars)
end#end integrate

