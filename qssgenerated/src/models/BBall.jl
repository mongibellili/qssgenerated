function f(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(-9.8, mulT(d[1], muladdT(100000.0, q[1], mulT(30.0, q[2], cache[4]), cache[3]), cache[2]), cache[1])
        #= none:1 =#
        return nothing
    end
end
