function changeVarNames_to_q_d(ex::Expr,stateVarName::Symbol)
  newEx=postwalk(ex) do a  # change rhs equ to contain q[] and d[] instead of user symbols
      if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
          a.args[1]!=stateVarName && a.args[1]!=:discrete && error("symbol $(a.args[1]) undefined")  # for now I do not allow parameters
          if a.args[1]==stateVarName 
              a.args[1]=:q 
           else
               a.args[1]=:d
           end
      end
      return a
  end
  newEx
end
function changeVarNames_to_x_d(ex::Expr,stateVarName::Symbol)######maybe x is better for zc...if thats the case use this inside
newEx=postwalk(ex) do a  # change rhs equ to contain q[] and d[] instead of user symbols
    if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
        a.args[1]!=stateVarName && a.args[1]!=:discrete && error("symbol $(a.args[1]) undefined")
        if a.args[1]==stateVarName 
            a.args[1]=:x 
         else
             a.args[1]=:d
         end
    end
    return a
end
newEx
end

function replace_parameters(ex::Expr,param::Dict{Symbol,Number})
  postwalk(ex) do a  # change rhs equ to contain actual values instead of param symbols
    if a isa Symbol   
      if haskey(param, a)
       a=param[a] 
      end # no need to throw error if something undef, RuntimeGeneratedFunction will: tested(LoadError: UndefVarError: parameter3 not defined)
    end
    return a
  end
end



function extractJac_from_equs(ex::Expr,T::Int,D::Int,usymbols::Vector{SymEngine.Basic},dsymbols::Vector{SymEngine.Basic},jac::Vector{Vector{SymEngine.Basic}},jacDiscrete::Vector{Vector{SymEngine.Basic}})
  m=postwalk(ex) do a   #after equs constructed, eliminate ref ; use new expr m since we still need z for equs below (below the caller of this)...postwalk
      if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
       a=Symbol((a.args[1]), (a.args[2]))  #a become An #needed for differentiation ...jacobians....
      end
      return a
  end
  jacArr = []
  jacDiscArr = []
  basi = convert(Basic, m)
  #extract jaco components
  for j = 1:T            
      coef = diff(basi, usymbols[j])
      push!(jacArr, coef)
  end
  for j = 1:D            
      coef = diff(basi, dsymbols[j])
      push!(jacDiscArr, coef)
     # @show coef
  end
  push!(jac, jacArr)
  push!(jacDiscrete, jacDiscArr)             
end

function twoInOneSimplecase(ex)#  it s easier and healthier to leave the big case alone in one prewalk (below) when returning the size of the distributed cahce
  ex=Expr(:call, :createT,ex)# case where rhs of eq is a number needed to change to taylor (in cahce) to be used for qss ||  #case where rhs is q[i]...reutrn cache that contains it
  cachexpr = Expr(:ref, :cache)   # add cache[1] to createT(ex,....)
  push!(cachexpr.args,1)
  push!( ex.args, cachexpr)    
  return ex
end

function twoInOne(ex)# name to be changed later....i call this funciton in the docs the function that does the transformation
 cachexpr_lengthtracker = Expr(:mongi)# an expr to track the number of distributed caches
 
  prewalk(ex) do x
    #############################minus sign#############################################
    if x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 3
            if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3
                push!(x.args, x.args[3])#  args[4]=c
                x.args[3] = x.args[2].args[3]# args[3]=b
                x.args[2] = x.args[2].args[2]# args[2]=a
                x.args[1] = :subsub
                push!(cachexpr_lengthtracker.args,:b) #b is anything ...not needed, just to fill vector tracker
                cachexpr = Expr(:ref, :cache) #prepare cache
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))#constrcut cache with index cache[1]
                push!(x.args, cachexpr)
            elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :+ && length(x.args[2].args) == 3
                push!(x.args, x.args[3])#  args[4]=c
                x.args[3] = x.args[2].args[3]# args[3]=b
                x.args[2] = x.args[2].args[2]# args[2]=a
                x.args[1] = :addsub # £ µ § ~....        
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                #cachexpr.args[2] = index[i]
                push!(x.args, cachexpr)
            elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
                push!(x.args, x.args[3])#  args[4]=c
                x.args[3] = x.args[2].args[3]# args[3]=b
                x.args[2] = x.args[2].args[2]# args[2]=a
                x.args[1] = :mulsub # £ µ § ~....
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr1 = Expr(:ref, :cache)
                push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr1)

            else
                x.args[1] = :subT  # symbol changed cuz avoid type taylor piracy    
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
            end
    elseif x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 2  && !(x.args[2] isa Number)#negate
            x.args[1] = :negateT  # symbol changed cuz avoid type taylor piracy
            push!(cachexpr_lengthtracker.args,:b)
            cachexpr = Expr(:ref, :cache)
            push!(cachexpr.args,length(cachexpr_lengthtracker.args))
            push!(x.args, cachexpr)
    ############################### plus sign#######################################
    elseif x isa Expr && x.head == :call && x.args[1] == :+ && length(x.args) == 3
          if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3        
                push!(x.args, x.args[3])#  args[4]=c
                x.args[3] = x.args[2].args[3]# args[3]=b
                x.args[2] = x.args[2].args[2]# args[2]=a
                x.args[1] = :subadd#:µ  # £  § ....
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                #cachexpr.args[2] = index[i]
                push!(x.args, cachexpr)
          elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
                push!(x.args, x.args[3])#  args[4]=c
                x.args[3] = x.args[2].args[3]# args[3]=b
                x.args[2] = x.args[2].args[2]# args[2]=a
                x.args[1] = :muladdT#:µ  # £  § ....
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr1 = Expr(:ref, :cache)
                push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr1)
          else
                x.args[1] = :addT
                push!(cachexpr_lengthtracker.args,:b)
                cachexpr = Expr(:ref, :cache)
                push!(cachexpr.args,length(cachexpr_lengthtracker.args))
                push!(x.args, cachexpr)
          end
    elseif x isa Expr && x.head == :call && x.args[1] == :+ && (4 <= length(x.args) <= 9) # special add :i stopped at 9 cuz by testing it was allocating anyway after 9
              x.args[1] = :addT
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
   ############################### multiply sign#######################################
      #never happens :#  elseif x isa Expr && x.head == :call && x.args[1]==:* && length(x.args)==3 to get addmul or submul 
    elseif x isa Expr && x.head == :call && (x.args[1] == :*) && (3 == length(x.args))
           x.args[1] = :mulT
           push!(cachexpr_lengthtracker.args,:b)
           cachexpr1 = Expr(:ref, :cache)
           push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
           push!(x.args, cachexpr1)
    elseif x isa Expr && x.head == :call && (x.args[1] == :*) && (4 <= length(x.args) <= 7)# i stopped at 7 cuz by testing it was allocating anyway after 7
            x.args[1] = :mulTT
            push!(cachexpr_lengthtracker.args,:b)
            cachexpr1 = Expr(:ref, :cache)
            push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
            push!(x.args, cachexpr1)
            push!(cachexpr_lengthtracker.args,:b)
            cachexpr2 = Expr(:ref, :cache)   #multiply needs two caches
            push!(cachexpr2.args,length(cachexpr_lengthtracker.args))
            push!(x.args, cachexpr2)
    elseif x isa Expr && x.head == :call && (x.args[1] == :/)
              x.args[1] = :divT   
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
    elseif x isa Expr && x.head == :call && (x.args[1] == :exp ||x.args[1] == :log ||x.args[1] == :sqrt ||x.args[1] == :abs )
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
                
    elseif x isa Expr && x.head == :call && (x.args[1] == :^ )
              x.args[1] = :powerT  # should be like above but for lets keep it now...in test qss.^ caused a problem...fix later
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
    elseif x isa Expr && x.head == :call && (x.args[1] == :cos ||x.args[1] == :sin ||x.args[1] == :tan||x.args[1] == :atan )
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
              #second cache
              push!(cachexpr_lengthtracker.args,:b)#index2
              cachexpr2 = Expr(:ref, :cache)   
              push!(cachexpr2.args,length(cachexpr_lengthtracker.args))#construct cache[2]
              push!(x.args, cachexpr2)
   
    elseif x isa Expr && x.head == :call && (x.args[1] == :acos ||x.args[1] == :asin)
              #did not change symbol here
              push!(cachexpr_lengthtracker.args,:b)
              cachexpr = Expr(:ref, :cache)
              push!(cachexpr.args,length(cachexpr_lengthtracker.args))
              push!(x.args, cachexpr)
              #second cache
              push!(cachexpr_lengthtracker.args,:b)#index2
              cachexpr2 = Expr(:ref, :cache)   
              push!(cachexpr2.args,length(cachexpr_lengthtracker.args))#construct cache[2]
              push!(x.args, cachexpr2)
              #third cache
              push!(cachexpr_lengthtracker.args,:b)#index3
              cachexpr2 = Expr(:ref, :cache)   
              push!(cachexpr2.args,length(cachexpr_lengthtracker.args))#construct cache[2]
              push!(x.args, cachexpr2)        
    elseif false #holder for "myviewing" for next symbol for other functions...
       
    end
    return x # this is the line that actually enables modification of the original expression (prewalk only)
  end#end prewalk

  #return ex ########################################************************
  ex.args[2]=length(cachexpr_lengthtracker.args)# return the number of caches distributed...later clean the max of all these

   return ex
# end #end if number else prewalk

end#end function

#this macro is to be deleted : created for testing
#= macro changeAST(ex)
  Base.remove_linenums!(ex)
  # dump(ex; maxdepth=18)
  twoInOne(ex)
 
# dump( ex.args[1])
#=   @show ex.args[1]# return  
return nothing =#
esc(ex.args[1])# return 
end =#