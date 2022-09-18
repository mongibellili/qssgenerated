using qssgenerated
using Plots;
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
   sol = QSS_Solve(odeprob,2.0,qss2())
  x=evaluateSol(sol,1,1.0)
  pointt=[1.0]
  pointx=[x]
  display(scatter(pointt, pointx,label="x1",line=(:dot, 6),color=:black))
   plotSol(sol)





  # plotSol(sol[1],sol[2])
  # plotSol(sol[1],sol[2][1])
  
   # plotSol(sol.savedTimes,sol.savedVars)
       #=  temp1 = []
        temp2 = []
        for i = 1:length(sol[2][1])
            push!(temp1, sol[2][1][i].coeffs[1])
            push!(temp2, sol[2][2][i].coeffs[1])
        end
        display(plot!(sol[1], temp1))
        display(plot!(sol[1], temp2))
       println("done")
        readline()  =#
   
end
test()
