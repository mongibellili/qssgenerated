using qssgenerated

include("D:\\qssgenerated\\qssgenerated\\qssv01\\src\\models\\classicProblem.jl") # where you saved the model
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
#=         du[1] = u[2]
        du[2] =-u[1]-u[2] =#
        #= du[1] = 0.01*u[2]
        du[2] =-100*u[1]-100*u[2]+2020 =#
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2
    end
  
   
    
    sol2=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,qss2()) # 
   # plotSol(sol)
    plotSol(sol2)
end
test()
