

using qssgenerated
#= using Plots;
gr(); =#
function test()
    odeprob = @NLodeProblem begin
        parameter1=30.0# cache can be dynamic....parameters take this feature
        parameter2=1e+5
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(parameter2*u[1]+parameter1*u[2])
        if u[1]>0   #5*discrte gave error
            discrete[1]=0.0   #discrete=0.0-->type Symbol has no field args...find to personalize error msg            
        else
            discrete[1]=1.0                                    
        end
    end
    sol= QSS_Solve(odeprob,5.0,qss2())
    #save_prob_to_model(odeprob) #this requires macro prob
    #sol= QSS_Solve_from_model(5.0,qss2()) # turn off macro prob
  #  plotSol(sol[1],sol[2][1])
       #=  temp1 = []
        temp2 = []
        for i = 1:length(sol[2][1])
            push!(temp1, sol[2][1][i].coeffs[1])
            push!(temp2, sol[2][2][i].coeffs[1])
        end
        display(plot!(sol[1], temp1))
        display(plot!(sol[1], temp2))
        println("done")
        readline() =#
end
test()
