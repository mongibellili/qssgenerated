

(addsub)(a, b, c)=(-)((+)(a,b),c)
(subsub)(a, b, c)=(-)((-)(a,b),c)
(subadd)(a, b, c)=(+)((-)(a,b),c)# this should not exist cuz addsub has 3 args + cache

  function addTfoldl(op, res, bs...)# op is only addt so remove later
      l = length(bs)
     
      #i =  0;  l == i && return a  # already checked
      l == 1 && return res 
      i =  2; l == i && return op(res, bs[1],bs[end])
      i =  3;  l == i && return op(res, bs[1], bs[2],bs[end])
      i =  4;  l == i && return op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[end])
      i =  5;  l == i && return op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end])
      i =  6;  l == i && return op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[end])
      i =  7;  l == i && return op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[6],bs[end])
      i =  8;  l == i && return op(op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[6],bs[end]),bs[7],bs[end])
      i =  9;y=op(op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[6],bs[end]),bs[7],bs[8],bs[end]);  l == i && return y
      #last i creates y and passes it to the next loop: note that from this point, op will allocate
      #  warn the users . (inserting -0 will avoid allocation lol!)
       for i in i:l-1 # -1 cuz last one is cache...but these lines are never reached cuz macro does not create mulT for args >7 in the first place
          y = op(y, bs[i],bs[end])
      end
      return y 
  end
    function addT(a, b, c,d, xs...)
      if length(xs)!=0# if ==0 case below where d==cache
     
        addTfoldl(addT, addT( addT(a,b,c,xs[end]) , d ,xs[end] ), xs...)
      end

    end

    function mulTfoldl(op, res, bs...) # op is only mulT so remove later
      l = length(bs)
    #  @show l
      #println("entered multfoldl")
      #i =  0;  l == i && return a  # already checked
      i =  2;  l == i && return res 
      i =  3;  l == i && return              op(res, bs[1],bs[end-1],bs[end])
      i =  4; l == i && return           op(op(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end])
      i =  5;  l == i && return        op(op(op(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end]),bs[3],bs[end-1],bs[end])


    end
    function mulT(a, b, c, xs...)
     # println("entered mult(a b c)")
      if length(xs)>1# length(xs)==0 will never happen, length(xs)==1 is the case far-below
        mulTfoldl(mulT, mulT( mulT(a,b,xs[end-1],xs[end]) , c,xs[end-1] ,xs[end] ), xs...)
      end
    end

#all these functions can be written in terms of each other without cost (I think). i will test later.

  # all methods should have new names . no type piracy!!!
  function addsub(a::Taylor0{T}, b::Taylor0{T},c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
   #@__dot__ cache.coeffs = addsub(a.coeffs, b.coeffs,c.coeffs)
   @__dot__ cache.coeffs=(-)((+)(a.coeffs,b.coeffs),c.coeffs)
    return cache
  end
  function addsub(a::Taylor0{T}, b::Taylor0{T},c::T,cache::Taylor0{T}) where {T<:Number}
            cache.coeffs.=b.coeffs  
            cache[0]=cache[0]-c
   @__dot__ cache.coeffs = (+)(cache.coeffs, a.coeffs)
    return cache
  end
  function addsub(a::T, b::Taylor0{T},c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
            cache.coeffs.=b.coeffs  
            cache[0]=a+ cache[0]
           # @show a
   @__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
    return cache
  end
  addsub(a::Taylor0{T}, b::T,c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=addsub(b, a ,c,cache)#addsub(b::T, a::Taylor0{T} ,c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
 
  function addsub(a::T, b::T,c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
    #cache.coeffs.=-c.coeffs #not tested
    #cache.coeffs.=0.0
    @__dot__ cache.coeffs = (-)(c.coeffs)
   # @__dot__ cache.coeffs=(-)(cache.coeffs, c.coeffs)#   ==-c
    cache[0]=a+b+ cache[0] 
    return cache
end



function addsub(c::Taylor0{T},a::T, b::T,cache::Taylor0{T}) where {T<:Number}
  cache.coeffs.=c.coeffs 
  cache[0]=a+ cache[0]-b
  return cache
end
 addsub(a::T, c::Taylor0{T},b::T,cache::Taylor0{T}) where {T<:Number}= addsub(c,a, b,cache) 
 
 function addsub(c::T,a::T, b::T,cache::Taylor0{T}) where {T<:Number}
  ##[[[[[[[[[!!!!!!!!!!!!!!!!!]]]]]]]]]cache[1]...should be emptied
  #cache.coeffs.=0.0   # since C is number this cant be an op in the middle (all in middle op return taylor)
                      # still the cache can be full from another completely diff equation
                      # emptying the cache does not affect c (in this case c is not cache)
                      #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
  cache[0]=a+c-b
  return cache
end

###########"negate###########
function negateT(a::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
  @__dot__ cache.coeffs = (-)(a.coeffs)
  return cache
end
#= function negateT(b::T,cache::Taylor0{T}) where {T<:Number} # no need since ast of -Number is direct
  cache.coeffs.=0.0
  cache[0]=-b
  return cache
end =#
#################################################subsub########################################################""


  function subsub(a::Taylor0{T}, b::Taylor0{T},c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
    
    @__dot__ cache.coeffs = subsub(a.coeffs, b.coeffs,c.coeffs)
    
    return cache
  end

  function subsub(a::Taylor0{T}, b::Taylor0{T},c::T,cache::Taylor0{T}) where {T<:Number}
    #@show a
   cache.coeffs.=b.coeffs  
   cache[0]=cache[0]+c     #(a-b-c=a-(b+c))
  @__dot__ cache.coeffs = (-)(a.coeffs, cache.coeffs)
    return cache
  end
  subsub(a::Taylor0{T}, b::T,c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=subsub(a, c,b,cache) 
   #------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------
 #


 #      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #--------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-------
  function subsub(a::T,b::Taylor0{T}, c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
   # @show c
# bug: this is the first time cache is used upfront..!!!cache not empty
#cache.coeffs.=0.0
@__dot__ cache.coeffs = (-)(b.coeffs) # like before since a is number modifying the cache is fine
  # @__dot__ cache.coeffs=(-)(cache.coeffs, b.coeffs)
   cache[0]=a+ cache[0]
  # @show a
@__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
return cache



   #=  cache.coeffs.=-b.coeffs #not tested
    cache[0]=a+cache[0]     #(a-b-c=a-(b+c))
   @__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
     return cache =#
   end


   function subsub(a::Taylor0{T}, b::T,c::T,cache::Taylor0{T}) where {T<:Number}
    cache.coeffs.=a.coeffs  # if first op then needed...else cache is already  a
                                #maybe later have firstOp for every op
    cache[0]=cache[0]-b-c     
     return cache
   end

   function subsub( a::T,b::Taylor0{T},c::T,cache::Taylor0{T}) where {T<:Number}
    #cache.coeffs.=-a.coeffs  #not tested
   # cache.coeffs.=0.0 #]]]]]]]]]]]]]]]]]]]]]]$$$$$$$$$$$$$$$$********************
    @__dot__ cache.coeffs = (-)(b.coeffs)
   # @__dot__ cache.coeffs=(-)(cache.coeffs, a.coeffs)
    cache[0]=cache[0]+a-c     
     return cache
   end
   subsub( b::T,c::T,a::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=subsub( b,a,c,cache) 


   function subsub( a::T,b::T,c::T,cache::Taylor0{T}) where {T<:Number}
    #cache.coeffs.=0.0 # a is number, ok to empty cache
     #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
    cache[0]=a-b-c    
     return cache
   end


   #################################subadd####################################""
  #= function subadd(a::Taylor0{T}, b::Taylor0{T},c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}

   @__dot__ cache.coeffs = subadd(a.coeffs, b.coeffs,c.coeffs)
    #return 0#$T(v, a.order)
    return cache
  end =#

  #subadd(a,b,c)::Taylor0{T} where {T<:Number}=addsub(a,c,b)   #not tested
  function subadd(a::P,b::Q,c::R,cache::Taylor0{T}) where {P,Q,R <:Union{Taylor0,Number},T<:Number}
    addsub(a,c,b,cache)  
end

  ##################################""subT################################
  function subT(a::Taylor0{T}, b::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
   @__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
    return cache
  end
  function subT(a::Taylor0{T}, b::T,cache::Taylor0{T}) where {T<:Number}
    #= for i=1:length(cache.coeffs)
        cache.coeffs[i]=a.coeffs[i] 
    end =#
    @__dot__ cache.coeffs = a.coeffs  # not needed for in middle ops
   cache[0]=cache[0]-b    
    #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     return cache
   end
   function subT(a::T, b::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
    #cache.coeffs.=-b.coeffs  #not tested
    #cache.coeffs.=0.0
    @__dot__ cache.coeffs = (-)(b.coeffs) # a is number ok to modify cache
   # @__dot__ cache.coeffs=(-)(cache.coeffs, b.coeffs)
   cache[0]=cache[0]+a    
    #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     return cache
   end
   function subT( a::T,b::T,cache::Taylor0{T}) where {T<:Number}
    #cache.coeffs.=0.0 # a is number ok to empty
     #i think emptying the cashe once before the equation is cheaper than emptying it 
                        #for every op
    cache[0]=a-b   
    #@show cache
     return cache
   end
function addT(a::Taylor0{T}, b::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
  
        @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
                return cache
end
function addT(a::Taylor0{T}, b::T,cache::Taylor0{T}) where {T<:Number}
  #@show a==cache
    cache.coeffs.=a.coeffs  #not needed for inmiddle ops
    cache[0]=cache[0]+b 
          return cache
end
addT(b::T,a::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=addT(a, b,cache) 

function addT( a::T,b::T,cache::Taylor0{T}) where {T<:Number}
  #cache.coeffs.=0.0 # a is number ok to empty
   #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
  cache[0]=a+b   
  #@show cache
   return cache
 end

 #add Three vars a b c
function addT(a::Taylor0{T}, b::Taylor0{T},c::Taylor0{T},cache::Taylor0{T}) where {T<:Number}
 # println("add abc a=",a)
 # println("add abc cache=",cache)
  
 # println("equal cache=",cache==a)

@__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs,c.coeffs)
  return cache
end   
function addT(a::Taylor0{T}, b::Taylor0{T},c::T,cache::Taylor0{T}) where {T<:Number}
  # badly thought....no need to put a in cache just cache=a+b and then cache0+=c
  cache.coeffs.=a.coeffs  #not needed for in middle ops
    cache[0]=cache[0]+c 
  @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
    return cache
  end   
  addT(a::Taylor0{T},c::T, b::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=addT(a, b,c,cache) 
  addT(c::T, a::Taylor0{T}, b::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=addT(a, b,c,cache)

  function addT(a::Taylor0{T}, b::T,c::T,cache::Taylor0{T}) where {T<:Number}
  
    cache.coeffs.=a.coeffs  #not needed for inmiddle ops
      cache[0]=cache[0]+c+b
      return cache
    end 
    addT( b::T,a::Taylor0{T},c::T,cache::Taylor0{T}) where {T<:Number}=addT(a, b,c,cache) 
    addT( b::T,c::T,a::Taylor0{T},cache::Taylor0{T}) where {T<:Number}=addT(a, b,c,cache) 

    function addT( a::T,b::T,c::T,cache::Taylor0{T}) where {T<:Number}
      #cache.coeffs.=0.0 #a is number ok to empty
       #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
      cache[0]=a+b+c    
       return cache
     end




####################mul uses one cache when the original mul has only two elemnts a * b
function mulT(a::Taylor0{T}, b::T,cache1::Taylor0{T}) where {T<:Number}
  fill!(cache1.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs * cache1.coeffs  ##fixed broadcast dimension mismatch
#maybe clear cache2 and see if it allocs
  return cache1
end
mulT(a::T,b::Taylor0{T}, cache1::Taylor0{T}) where {T<:Number} = mulT(b , a,cache1)

function mulT(a::Taylor0{T}, b::Taylor0{T},cache1::Taylor0{T}) where {T<:Number}
 for k in eachindex(a)
    #println("inloop a= ",a)
   # 
      @inbounds cache1[k] = a[0] * b[k]
      @inbounds for i = 1:k
        cache1[k] += a[i] * b[k-i]
      end
     # println("inloop cache= ",cache)
  end

 #@__dot__ cache1.coeffs = cache2.coeffs
#maybe clear cache2 and see if it allocs
  return cache1
end
function mulT(a::T, b::T,cache1::Taylor0{T}) where {T<:Number}
  #cache1 needs to be clean: a clear here will alloc...this is the only "mul" that wants a clean cache
  #sometimes the cache can be dirty at only first position: it will not throw an assert error!!!
 # cache[1].coeffs.=0.0 #it causes two allocs for mul(cst,cst,a,b)
  cache1[0]=a*b  
       return cache1
end
#########################mul uses two caches when the op has many terms###########################

function mulT(a::Taylor0{T}, b::T,cache1::Taylor0{T},cache2::Taylor0{T}) where {T<:Number}
  fill!(cache2.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs * cache2.coeffs  ##fixed broadcast dimension mismatch
#maybe clear cache2 and see if it allocs
  return cache1
end
mulT(a::T,b::Taylor0{T}, cache1::Taylor0{T},cache2::Taylor0{T}) where {T<:Number} = mulT(b , a,cache1,cache2)

function mulT(a::Taylor0{T}, b::Taylor0{T},cache1::Taylor0{T},cache2::Taylor0{T}) where {T<:Number}
 for k in eachindex(a)
    #println("inloop a= ",a)
   # 
      @inbounds cache2[k] = a[0] * b[k]
      @inbounds for i = 1:k
        cache2[k] += a[i] * b[k-i]
      end
     # println("inloop cache= ",cache)
  end

 @__dot__ cache1.coeffs = cache2.coeffs
#maybe clear cache2 and see if it allocs
  return cache1
end
function mulT(a::T, b::T,cache1::Taylor0{T},cache2::Taylor0{T}) where {T<:Number}
  #cache1 needs to be clean: a clear here will alloc...this is the only "mul" that wants a clean cache
  #sometimes the cache can be dirty at only first position: it will not throw an assert error!!!
 # cache[1].coeffs.=0.0 #it causes two allocs for mul(cst,cst,a,b)
  cache1[0]=a*b  
       return cache1
end


#########################################muladdT and subadd not test  : added T to muladd cuz there did not want to import muladd to extend...maybe later
function muladdT(a::Taylor0{T}, b::Taylor0{T}, c::Taylor0{T},cache1::Taylor0{T}) where {T<:Number}
  for k in eachindex(a)
       #@inbounds 
       cache1[k] = a[0] * b[k] + c[k]
       @inbounds for i = 1:k
         cache1[k] += a[i] * b[k-i]
       end
   end
  @__dot__ cache1.coeffs = cache2.coeffs
   return cache1
 end
#the rest of muladd and subadd to be implemented later
function muladdT(a::P,b::Q,c::R,cache1::Taylor0{T}) where {P,Q,R <:Union{Taylor0,Number},T<:Number}
      addT(mulT(a, b,cache1),c,cache1) # improvement of added performance not tested: there is one extra step of 
                                       #puting cache in itself
 end


 function mulsub(a::P,b::Q,c::R,cache1::Taylor0{T}) where {P,Q,R <:Union{Taylor0,Number},T<:Number}
  subT(mulT(a, b,cache1),c,cache1)
end

function divT(a::Taylor0{T}, b::T,cache1::Taylor0{T}) where {T<:Number}
  fill!(cache1.coeffs, b)
  #= @inbounds aux = a.coeffs[1] / b
  v = Array{typeof(aux)}(undef, length(a.coeffs)) =#
  @__dot__ cache1.coeffs = a.coeffs / cache1.coeffs
  return cache1
end

#divT(a::Taylor0{T}, b::Taylor0{S}) where {T<:Number,S<:Number} = divT(promote(a,b)...)

function divT(a::Taylor0{T}, b::Taylor0{T},cache1::Taylor0{T}) where {T<:Number}
  iszero(a) && !iszero(b) && return cache1 # this op has its own clean cache2
  #= if a.order != b.order
      a, b = fixorder(a, b)
  end =#
  # order and coefficient of first factorized term
  ordfact, cdivfact = divfactorization(a, b)
  cache1[0] = cdivfact
  #c = Taylor0(cdivfact, a.order-ordfact)
  #for ord in eachindex(c)
     # div!(c, a, b, ord) # updates c[ord]
  for k=1:a.order-ordfact
      imin = max(0, k+ordfact-b.order)
      @inbounds cache1[k] = cache1[imin] * b[k+ordfact-imin]
      @inbounds for i = imin+1:k-1
        cache1[k] += cache1[i] * b[k+ordfact-i]
      end
      if k+ordfact ≤ b.order
          @inbounds cache1[k] = (a[k+ordfact]-cache1[k]) / b[ordfact]
      else
          @inbounds cache1[k] = - cache1[k] / b[ordfact]
      end
  end
 # @__dot__ cache1.coeffs = cache2.coeffs
#maybe clear cache2 and see if it allocs
  return cache1
end


function clearCache(cache::Vector{Taylor0{Float64}},size::Int) #where {T<:Number}
  for i=1:size
    cache[i].coeffs.=0.0
  end
end