function f(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(negateT(q[1], cache[2]), q[2], cache[1])
        #= none:1 =#
        return nothing
    end
end
function twoVarSys1(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(negateT(q[1], cache[2]), q[2], cache[1])
        #= none:1 =#
        return nothing
    end
end
function twoVarSys12(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(negateT(q[1], cache[2]), q[2], cache[1])
        #= none:1 =#
        return nothing
    end
end
