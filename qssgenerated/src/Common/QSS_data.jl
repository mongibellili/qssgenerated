#hold helper datastructures needed for simulation, can be seen as the model in the qss architecture (model-integrator-quantizer)
struct QSS_data{T,Z}
    quantum :: Vector{Float64} 
    x :: Vector{Taylor0{Float64}}
    q :: Vector{Taylor0{Float64}}
    tx ::  MVector{T,Float64} 
    tq :: MVector{T,Float64} 
    nextStateTime :: MVector{T,Float64}    
    nextInputTime :: Vector{Float64}  # later change to Mvector
    nextEventTime :: MVector{Z,Float64}  
    t::Taylor0{Float64}# mute taylor var to be used with math functions to represent time
    integratorCache::Taylor0{Float64}
    order::Int
    savedVars :: Vector{Array{Taylor0{Float64}}} #has to be vector (not SA) cuz to be resized in integrator
    savedTimes :: Vector{Float64}  
    taylorOpsCache::Vector{Taylor0{Float64}}

    finalTime:: Float64   
    savetimeincrement::Float64 
    initialTime :: Float64    
    dQmin ::Float64    
    dQrel ::Float64  
end
function saveat(savetimeincrement::Float64) # it s better than the user entre a number...fool-proof
    savetimeincrement
end 
qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss1()=Val(4)
liqss2()=Val(5)

function getOrderfromSolverMethod(::Val{V}) where {V}  # @generated and inline did not enhance performance  
    if V==1 || V==2 || V==3
        return V
    else
        return V-3
    end
end


function save_prob_to_model(prob::NLODEProblem{T,D,Z,Y},path::String) where {T,D,Z,Y}
    open(path, "a") do io        #if no model name is given, default to name f
         println(io,string(prob.eqs))  
     end
end
function save_prob_to_model(prob::NLODEProblem{T,D,Z,Y},path::String,modelName::String) where {T,D,Z,Y}#user wants to change the name of the model
    newFunName=Symbol(modelName)
    isdefined(Main,newFunName) && error("model exists!")# if the user included the path, if not it will just append a second function with same name
    def=splitdef(prob.eqs)
    def[:name] = newFunName 
    newFun=combinedef(def)
    open(path, "a") do io        
         println(io,string(newFun))  
     end
end

function QSS_Solve_from_model(f::Function,prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},savetimeincrement::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64) where {T,D,Z,Y,V}
    QSS_Unique_Solve(f,prob,finalTime,Val(V),savetimeincrement,initialTime,dQmin,dQrel)
end

function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},savetimeincrement::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64) where {T,D,Z,Y,V}
     f=@RuntimeGeneratedFunction(prob.eqs)
     QSS_Unique_Solve(f,prob,finalTime,Val(V),savetimeincrement,initialTime,dQmin,dQrel)
 end
 
function QSS_Unique_Solve(f::Function,prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},savetimeincrement::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64) where {T,D,Z,Y,V}    
     order=getOrderfromSolverMethod(Val(V))# this is to make a difference b/w qss1 and liqss1 but if performance then have qss_solve for liqss?
     quantum =  zeros(T)
     x = Vector{Taylor0{Float64}}(undef, T)
     q = Vector{Taylor0{Float64}}(undef, T)
     nextStateTime = @MVector zeros(T)
     nextInputTime =  zeros(T)
     tx = @MVector zeros(T)
     tq = @MVector zeros(T)
     nextEventTime=@MVector zeros(Z)
     t = Taylor0(zeros(order + 1), order)
     t[1]=1.0
     t[0]=initialTime
     integratorCache=Taylor0(zeros(order+1),order) #for integratestate only
     sizehint=floor(Int64, 1.0+(finalTime/savetimeincrement)) ## this number can be found from ft and saveat (ft/saveat=5/0.1=50)+1 and maybe *2/3 factor (time can jump past a saving time) of user or default saveat if not provided
     savedVars = Vector{Array{Taylor0{Float64}}}(undef, T)# has to be vector (not SA) cuz to be resized in integrator
     for i = 1:T
         temparr=Array{Taylor0{Float64}}(undef, sizehint)
         for j=1:sizehint          
             temparr[j]=Taylor0(zeros(order+1),order) 
         end
         savedVars[i]=temparr
         x[i]=Taylor0(zeros(order + 1), order) 
         x[i][0]= prob.initConditions[i]
         q[i]=Taylor0(zeros(order+1), order)#q normally 1order lower than x but since we want f(q) to  be a taylor that holds all info (1,2,3), lets have q of same order and not update last coeff        
         tx[i] = initialTime
         tq[i] = initialTime
     end
     savedTimes = zeros(sizehint) #estimate like sizehint...later stiffness hint...to be multiplied by a stiffnessfactor
     savedTimes[1]=initialTime
     cacheSize=prob.cacheSize
     taylorOpsCache=Array{Taylor0{Float64},1}()
     for i=1:cacheSize
       push!(taylorOpsCache,Taylor0(zeros(order+1),order))
     end
     qssdata= QSS_data(quantum,x,q,tx,tq,nextStateTime,nextInputTime ,nextEventTime , t, integratorCache,order,savedVars,savedTimes,taylorOpsCache,finalTime,savetimeincrement, initialTime,dQmin,dQrel)
    
  #=    zcf = Vector{Function}()
     for i = 1:length(prob.zceqs)# later change to numberZC
       push!(zcf, @RuntimeGeneratedFunction(prob.zceqs[i].args[2])) #args[2] cuz there is extra stuff
     end =#


    if V==1 || V ==2 || V ==3
     QSS_integrate(Val(V),qssdata,prob,f)   # solver=Val(\d) needed to specialize the integrator...and f has to be passed alone, it cannot be burried in qss_data...x2 performance
    else
     LiQSS_integrate(Val(V-3),qssdata,prob,f)
    end
     #return nothing be careful to add this
 end
 
 function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},savetimeincrement::Float64) where {T,D,Z,Y,V}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3
    QSS_Solve(prob,finalTime,Val(V),savetimeincrement,initialTime,dQmin,dQrel)
end

function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V}) where {T,D,Z,Y,V}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
    QSS_Solve(prob,finalTime,Val(V),savetimeincrement,initialTime,dQmin,dQrel)
end
function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64) where {T,D,Z,Y}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
    QSS_Solve(prob,finalTime,Val(1),savetimeincrement,initialTime,dQmin,dQrel)
end
function QSS_Solve(prob::NLODEProblem{T,D,Z,Y}) where {T,D,Z,Y}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1;finalTime=5.0
    QSS_Solve(prob,finalTime,Val(1),savetimeincrement,initialTime,dQmin,dQrel)
end

function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},savetimeincrement::Float64) where {T,D,Z,Y,V}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3
    QSS_Solve_from_model(f,m,finalTime,Val(V),savetimeincrement,initialTime,dQmin,dQrel)
    end

function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V}) where {T,D,Z,Y,V}
initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
QSS_Solve_from_model(f,m,finalTime,Val(V),savetimeincrement,initialTime,dQmin,dQrel)
end
function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y},finalTime::Float64) where {T,D,Z,Y}
initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
QSS_Solve_from_model(f,m,finalTime,Val(1),savetimeincrement,initialTime,dQmin,dQrel)
end
function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y}) where {T,D,Z,Y}
initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1;finalTime=5.0
QSS_Solve_from_model(f,m,finalTime,Val(1),savetimeincrement,initialTime,dQmin,dQrel)
end



