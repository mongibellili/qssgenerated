

using qssgenerated

function test()
    odeprob = @NLodeProblem begin
        cacheT=10
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(1e+5*u[1]+30.0*u[2])
        if u[1]>0   #5*discrte gave error
            discrete[1]=0.0               
        else
            discrete[1]=1.0                                    
        end
    end
    save_prob_to_model(odeprob,"/home/unknown/qssgenerated/qssgenerated/src/models/BBall.jl")        
end
test()
