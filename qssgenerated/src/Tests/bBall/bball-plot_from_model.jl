

using qssgenerated
using Plots;
include("/home/unknown/qssgenerated/qssgenerated/src/models/BBall.jl")# the location where you saved
function test()
    odeprob = @NLodeProblem begin
        k=1e+5
        g=9.8
        b=30.0
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-g-(discrete[1])*(k*u[1]+b*u[2])
        if u[1]>0   
            discrete[1]=0.0               
        else
            discrete[1]=1.0                                    
        end
    end
    sol=QSS_Solve_from_model(f,odeprob,8.0,qss2())
    x=evaluateSol(sol,1,1.0)
    pointt=[1.0]
    pointx=[x]
    display(scatter(pointt, pointx,label="x1",line=(:dot, 6),color=:black))
    plotSol(sol)
end
test()
