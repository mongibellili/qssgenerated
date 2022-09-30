using OrdinaryDiffEq
using qssgenerated
using BenchmarkTools
#using Plots
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]
        du[2] =-u[1]-u[2]
       #=  du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#
    
    end
    tspan = (0.0,5)
    u0 = [1.0,2.4]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
    #sol = solve(prob,Rosenbrock23(),abstol = 1e-6, reltol = 1e-3)
    
   # display(sol)
    #display(plot!(sol,line=(:dot, 4)))
   #=  display(plot!(sol,xlims=(0.0,0.5),ylims=(0.9,1.1)))
    println("done")
    readline() =#
end
@btime odeDiffEquPackage()  




#= function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        #=         du[1] = u[2]
        du[2] =-u[1]-u[2] =#
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0
        #= du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
    end
   sol = QSS_Solve(odeprob,2.0,liqss2())
   plotSol(sol)
 
   
end
test() =#
#odeDiffEquPackage() 