using qssgenerated
using BenchmarkTools
include("D:\\qssv01\\qssv0.1\\qssv01\\src\\models\\classicProblem.jl")


function test()
     odeprob = @NLodeProblem begin
          u = [1.0, 0.0]
          discrete = [0.0]
          du[1] = u[2]
          du[2] =-u[1]-u[2]
          #= du[1] = 0.01*u[2]
          du[2] =-100*u[1]-100*u[2]+2020 =#
          #= du[1] = -u[1]-u[2]+0.2
          du[2] =u[1]-u[2]+1.2 =#
      end
<<<<<<< HEAD
     sol=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,qss2())
=======
     sol=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss1())
>>>>>>> f6146dcfae368a3c045fb2f48b5714650e0e31dc
end
@btime test()
