
using BenchmarkTools
using qssgenerated
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
         #=  du[1] = 0.01*u[2]
          du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#
          #= du[1] = -u[1]-u[2]+0.2
          du[2] =u[1]-u[2]+1.2 =#
    end
   sol = QSS_Solve(odeprob,5.0,liqss3())
end
@btime test()
