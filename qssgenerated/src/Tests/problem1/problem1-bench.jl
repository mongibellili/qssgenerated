
using BenchmarkTools
using qssgenerated
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
   sol = QSS_Solve(odeprob,5.0,qss2())
end
@btime test()
