using qssgenerated
using BenchmarkTools
include("/home/unknown/qssgenerated/qssgenerated/src/models/classicProblem.jl") 


function test()
     odeprob = @NLodeProblem begin
          u = [1.0, 0.0]
          discrete = [0.0]
          du[1] = u[2]
          du[2] =-u[1]-u[2]
      end
     sol=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss1())
end
@btime test()
