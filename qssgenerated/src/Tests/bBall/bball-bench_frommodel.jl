
using BenchmarkTools
using qssgenerated
include("/home/unknown/qssgenerated/qssgenerated/src/models/classicProblem.jl") 
function test()
    odeprob = @NLodeProblem begin
        cacheT=10
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(1e+5*u[1]+30.0*u[2])
        if u[1]>0   #5*discrte gave error
            discrete[1]=0.0               
        else
            discrete[1]=1.0  
                                  
        end     
    end
    sol = QSS_Solve_from_model(f,odeprob,5.0,qss2())
    #save_prob_to_model(odeprob) #this requires macro prob
    #sol= QSS_Solve_from_model(5.0,qss2()) # turn off macro prob
end
@btime test()





































#=  using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]#*cos(t) #*30*exp(t)
        du[2] = -u[1] - u[2]#+sqrt(t) # +1#+10-t*t+t +cos(t)
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)

    end
    tspan = (0.0,5)
    u0 = [1.0,2.4]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
odeDiffEquPackage() =#



#@btime modqss.qss_Integrate(initCond,jacobian,order)


        #output=
#=         jacobian= StaticArrays.SVector{2, SymEngine.Basic}[[0, 1], [-100000.0*d1, -30.0*d1]]
        discJac= StaticArrays.SVector{1, SymEngine.Basic}[[0], [-100000.0*q1 - 30.0*q2]]
        zc_jac= StaticArrays.SVector{2, SymEngine.Basic}[[1, 0]]
        ZC_discJac= StaticArrays.SVector{1, SymEngine.Basic}[[0]]
        evHandlr= qss.EventHandlerStruct[qss.EventHandlerStruct{2, 1}(1, [NaN, NaN], [0.0]), qss.EventHandlerStruct{2, 1}(2, [NaN, NaN], [1.0])]
        SD= StaticArrays.SVector{2, Int64}[[0, 2],[1, 2] ]
        SZ= StaticArrays.SVector{1, Int64}[[1], [0]]
        HZ= StaticArrays.SVector{1, Int64}[[0], [0]]
        HD= StaticArrays.SVector{2, Int64}[[0, 2], [0, 2]] =#
