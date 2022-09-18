using qssgenerated

include("/home/unknown/qssgenerated/qssgenerated/src/models/classicProblem.jl") # where you saved the model
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 2.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
  
    sol=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,qss3(),saveat(0.1),0.0,1e-6,1e-3) # later if to change f to meaningfull name and maybe save all models in one file, user has to enter model name in the macro so that f be name upfront
    plotSol(sol)
end
test()
