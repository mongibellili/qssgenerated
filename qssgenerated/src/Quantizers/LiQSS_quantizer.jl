function updateQ(::Val{1},i::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::Vector{Array{Float64}},u::Vector{Float64})#####
q[i][0]=x[i][0]
q_=q[i][0]
dq=0.0
if x[i][1]>0
    dx=u[i]+(q_+quantum[i])*a[i][i]
    if dx>=0
        dq=quantum[i]
    else
        dq=(-u[i]/a[i][i])-q_
    end
else
    dx=u[i]+(q_-quantum[i])*a[i][i]
    if dx<=0
        dq=-quantum[i]
    else
        dq=(-u[i]/a[i][i])-q_
    end
end
if dq>2*quantum[i]
    dq=2*quantum[i]
end
if dq<-2*quantum[i]
    dq=-2*quantum[i]
end
q[i][0]=q_+dq
end
function Liqss_reComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    dt=0.0
    q=qv[i][0]
    x=xv[i][0]
    if xv[i][1] !=0.0
        dt=(q-x)/xv[i][1]
        if dt>0.0
            nextTime[i]=currentTime+dt
        else
            if xv[i][1]>0.0  
                nextTime[i]=currentTime+(q-x+2*quantum[i])/xv[i][1]
            else
                nextTime[i]=currentTime+(q-x-2*quantum[i])/xv[i][1]
            end
            if nextTime[i] < currentTime
                nextTime[i]=currentTime+1e-6
            end
        end
    else
        nextTime[i]=Inf
    end
end



function updateLinearApprox(::Val{1},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::Vector{Array{Float64}},u::Vector{Float64},qOld::Float64,derxOld::Float64)
    diffQ=q[i][0]-qOld
    if diffQ != 0.0
        a[i][i]=(x[i][1]-derxOld)/diffQ
    else
        a[i][i]=0.0
    end
    u[i]=x[i][1]-a[i][i]*q[i][0]
end