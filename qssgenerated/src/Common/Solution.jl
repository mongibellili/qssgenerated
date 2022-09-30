
struct Sol
    savedTimes::Vector{Float64} 
    savedVars::Vector{Array{Taylor0{Float64}}}
end
function plotSol(sol::Sol)
    numPoints=length(sol.savedTimes)
    numVars=length(sol.savedVars)
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, sol.savedVars[k][i].coeffs[1])
      end
      display(plot!(sol.savedTimes, temp,label="x$k"))    
    end
      println("press enter to exit")
      readline()
end
function plotSol(savedTimes::Vector{Float64} ,  savedVars::Vector{Array{Taylor0{Float64}}})
    numPoints=length(savedTimes)
    numVars=length(savedVars)
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, savedVars[k][i].coeffs[1])
      end
      display(plot!(savedTimes, temp,label="x$k"))
    end
      println("press enter to exit")
      readline()
end
function plotSol(savedTimes::Vector{Float64} , savedVar::Vector{Taylor0{Float64}})
      numPoints=length(savedTimes)
      temp = []
      for i = 1:numPoints
          push!(temp, savedVar[i].coeffs[1])
      end
      display(plot!(savedTimes, temp))
      println("press enter to exit")
      readline()
end
function evaluateSol(sol::Sol,index::Int,t::Float64)
        for i=1:length(sol[1])#savedTimes
            if sol[1][i]>t # i is closest lower point
                return sol[2][index][i-1](t-sol[1][i-1])#taylor evaluation after small elapsed with the point before (i-1)
            end
        end
end
(sol::Sol)(index::Int,t::Float64) = evaluateSol(p,index,t)

  
function getindex(s::Sol, i::Int64)
      if i==1
         return s.savedTimes
      elseif i==2
         return s.savedVars
      else
         error("sol has 2 attributes: time and states")
      end
end