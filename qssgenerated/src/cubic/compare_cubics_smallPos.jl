using LinearAlgebra
function linear(a::F, b::F) where F <: AbstractFloat
  if a == 0.0
    res = Inf
  else
    res = -b / a
    if res<=0
      res=Inf
    end
  end
  return res
end
function allreallinear(a::F, b::F) where F <: AbstractFloat
  if a == 0.0
    res = Inf
  else
    res = -b / a
  end
  return res
end
function allrealquadratic(a::F, b::F, c::F) where F <: AbstractFloat
  if a == 0.0
    res1=Inf
    res2 = allreallinear(b,c)
  else
    if c == 0.0
      res1=0.0
      res2 = allreallinear(a,b)
    else
      if b == 0.0
        r = -c / a
        if r <= 0.0
          res1=Inf
          res2=Inf
        else
          res1 = sqrt(r)
          res1 = -res1
        end
      else
        Δ = 1.0 - 4c*a / (b*b)
        if Δ < 0.0
          res1=Inf
          res2=Inf
        else
          q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
          res1 = q / a
         
          res2=c / q
          
        end
      end
    end
  end
  return res1,res2
end
function quadratic(a::F, b::F, c::F) where F <: AbstractFloat
  if a == 0.0
    res = linear(b,c)
  else
    if c == 0.0
      res = linear(a,b)
    else
      if b == 0.0
        r = -c / a
        if r <= 0.0
          res=Inf
        else
          res = sqrt(r)
        end
      else
        Δ = 1.0 - 4c*a / (b*b)
        if Δ < 0.0
          res=Inf
        else
          q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
          res = q / a
          if res<=0
            res=Inf
          end
          res2=c / q
          if  res2>0 && res2<res
            res=res2
          end
        end
      end
    end
  end
  return res
end
function cubic1(a::F, b::F, c::F, d::F) where F <: AbstractFloat
            if a == 0.0
                  res = quadratic(b,c,d)
            else
                      if d == 0.0
                            res = quadratic(a,b,c)
                      else
                                    A = b/a
                                    B = c/a
                                    C = d/a
                                    Q = (A^2-3B)/9
                                    R = (2*A^3-9A*B+27C)/54
                                    S = -A/3
                                    if R^2 < Q^3
                                            P = -2*sqrt(Q)
                                            ϕ = acos(R/sqrt(Q^3))
                                            res = P*cos(ϕ/3)+S
                                            @show res
                                            if res<=0
                                              res=Inf
                                            end
                                            res2 = P*cos((ϕ+2π)/3)+S
                                            #@show res2
                                            if res2>0 && res2 <res
                                              res=res2
                                            end
                                            res3 = P*cos((ϕ-2π)/3)+S
                                            #@show res3
                                            if res3>0 && res3 <res
                                              res=res3
                                            end
                                    else
                                            T = -sign(R)*cbrt(abs(R)+sqrt(R^2-Q^3))
                                            U = 0.0
                                            if T != 0.0
                                              U = Q/T
                                            end
                                            res= 0.5*(T+U)
                                            println("1root= ",res)
                                            if res<=0
                                              res=Inf
                                            end
                                    end
                      end
            end
  return res
end
function cubic2(a::F, b::F, c::F, d::F) where F <: AbstractFloat
  if a == 0.0
    res = quadratic(b,c,d)
  else
    if d == 0.0
      res = quadratic(a,b,c)
    else
      A = b/a
      B = c/a
      C = d/a
      Q = (A^2-3B)/9
      R = (2*A^3-9A*B+27C)/54
      p=-3*Q
      q=2*R
      S = -A/3
      if 4p*p*p+27q*q>0
        m=sqrt((q*q/4)+p*p*p/27)
        c=cbrt(-q/2+m)
        res=c-p/(3c)+S
        #if <0 Inf
      else
        res=Inf
      end

    end
  end
  return res
end
function cubic3(a::F, b::F, c::F, d::F) where F <: AbstractFloat
  if a == 0.0
    res = quadratic(b,c,d)
  else
    if d == 0.0
      res = quadratic(a,b,c)
    else
      A = b/a
      B = c/a
      C = d/a
      Q = (A^2-3B)/9
      R = (2*A^3-9A*B+27C)/54
      p=-3*Q
      q=2*R
      S = -A/3
      if 4p*p*p+27q*q>0
        m=sqrt((q*q/4)+p*p*p/27)
        res=cbrt(-q/2+m)+cbrt(-q/2-m)+S
      else
        res=Inf
      end
    end
  end
  return res
end
function cubic4(a::F, b::F, c::F, d::F) where F <: AbstractFloat#Cardano’s method:
  if a == 0.0
    res = quadratic(b,c,d)
  else
      if d == 0.0
              res = quadratic(a,b,c)
      else
              A = b/a
              B = c/a
              C = d/a
              Q = (A^2-3B)/9
              R = (2*A^3-9A*B+27C)/54
              p=-3*Q
              q=2*R
              S = -A/3
              if 4p*p*p+27q*q>0
                    resVect = allrealquadratic(1.0,q,-(p^3)/27)
                    z=resVect[1]
                  #  @show resVect
                    α=cbrt(resVect[1])
                    β=-p/(3α)
                    res=α+β+S 
                    if res<=0
                      res=Inf
                    end
                    α=cbrt(resVect[2])
                    β=-p/(3α)
                    res2=α+β+S
                    if res2 < res && res2 >0
                      res=res2
                    end 
              else
                  res=Inf        
              end
      end
  end
  return res#,res2#,res3
end
function cubic5(a::Float64, b::Float64, c::Float64, d::Float64)::Float64 #where F <: AbstractFloat
  _a = 1.0 / a
  b, c, d = b * _a, c * _a, d * _a
  m = b < c ? b : c
  m = d < m ? d : m
  m > eps(Float64) && return Inf#typemax(Float64) # Cauchy bound
  _3 = 1.0 / 3
  _9 = 1.0 / 9
  SQ3 = sqrt(3.0)
  xₙ = -b * _3
  b²_9 = b * b * _9
  yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)   #eq to 2R
  δ² = muladd(-_3, c, b²_9)                  #eq to Q
  h² = 4δ² * δ² * δ²
  Δ = muladd(yₙ, yₙ, -h²)
  if Δ > 4eps(Float64) # one real root and two complex roots
  p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
  q = δ² / p
  z = xₙ + p + q
  z > -eps(Float64) ? z : Inf#typemax(Float64)
  elseif Δ < -4eps(Float64) # three real roots
  θ = abs(yₙ) < eps(Float64) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
  δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
  z₁ = 2δ * cos(θ)
  z₂ = muladd(-0.5, z₁, xₙ)
  z₃ = SQ3 * δ * sin(θ)
  x₁ = xₙ + z₁
  x₂ = z₂ + z₃
  x₃ = z₂ - z₃
  x = x₁ > -eps(Float64) ? x₁ : Inf# typemax(F)
  x = x₂ > -eps(Float64) && x₂ < x ? x₂ : x
  x₃ > -eps(Float64) && x₃ < x ? x₃ : x
  else # double or triple real roots
  δ = cbrt(0.5yₙ)
  x₁ = xₙ + δ
  x₂ = xₙ - 2δ
  x = x₁ > -eps(Float64) ? x₁ : Inf#typemax(F)
  x₂ > -eps(Float64) && x₂ < x ? x₂ : x
  end
end
using BenchmarkTools
function closeEnough(sum1::Float64,sum2::Float64) 
    return abs(sum1-sum2)>1e-12
end
#= @show cubic1(-1.33,2.0,1.0,6.0)#1.007151552710968
@show cubic2(-1.33,2.0,1.0,6.0) #2.515556238254016
@show cubic3(-1.33,2.0,1.0,6.0) #2.515556238254017
@show cubic4(-1.33,2.0,1.0,6.0) # 2.5155562382540193
@show cubic5(-1.33,2.0,1.0,6.0)#2.515556238254016 =#

#= @show cubic1(1.33,2.0,1.0,6.0)#-0.8186606842427812
@show cubic2(1.33,2.0,1.0,6.0) #Inf
@show cubic3(1.33,2.0,1.0,6.0) #2.515556238254017
@show cubic4(1.33,2.0,1.0,6.0) # 2.5155562382540193
@show cubic5(1.33,2.0,1.0,6.0)#2.515556238254016 =#



a=-1.33
b=2.0
c=1.0
d=6.0
  
#(-1.33,2.0,1.0,6.0)
#@btime cubic5(1.33,12.65,1.0,-6.0)#   0.6308983740365341 #63.615 ns (0 allocations: 0 bytes)
@show cubic5(a,b,c,d)
@btime cubic5(a,b,c,d)#               0.6308983740365341 #87.923 ns (1 allocation: 16 bytes)#69.736 ns (0 allocations: 0 bytes) if  change f to float64
#@show cubic5(-1.33,2.0,1.0,6.0)      #2.515556238254016  #19.098 ns (0 allocations: 0 bytes)