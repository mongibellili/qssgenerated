using LinearAlgebra

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

function cubic1(a::F, b::F, c::F, d::F) where F <: AbstractFloat
            if a == 0.0
                  res1,res2 = allrealquadratic(b,c,d)
                  res3=Inf
            else
                if d == 0.0
                  res1,res2 = allrealquadratic(a,b,c)
                  res3=0.0
                else
                    A = b/a
                    B = c/a
                    C = d/a
                    Q = (A^2-3B)/9
                    R = (2*A^3-9A*B+27C)/54
                    S = -A/3
                    if R^2 < Q^3   #3 real
                            P = -2*sqrt(Q)
                            ϕ = acos(R/sqrt(Q^3))
                            res1 = P*cos(ϕ/3)+S
                            
                            res2 = P*cos((ϕ+2π)/3)+S
                            
                            res3 = P*cos((ϕ-2π)/3)+S
                            
                    else   #2Imag 1 real
                            T = -sign(R)*cbrt(abs(R)+sqrt(R^2-Q^3))
                            U = 0.0
                            if T != 0.0
                              U = Q/T
                            end
                            res1= 0.5*(T+U)
                            res2=Inf
                            res3=Inf                                           
                    end
                end
            end
  return res1,res2,res3
end
function cubic2(a::F, b::F, c::F, d::F) where F <: AbstractFloat
  if a == 0.0
    res1,res2 = allrealquadratic(b,c,d)
    res3=Inf
  else
    if d == 0.0
      res1,res2 = allrealquadratic(a,b,c)
      res3=0.0
    else
      A = b/a
      B = c/a
      C = d/a
      Q = (A^2-3B)/9
      R = (2*A^3-9A*B+27C)/54
      p=-3*Q
      q=2*R
      S = -A/3
      if 4p*p*p+27q*q>0  #2Imag 1 real
        m=sqrt((q*q/4)+p*p*p/27)
        c=cbrt(-q/2+m)
        res1=c-p/(3c)+S
        res2=Inf
        res3=Inf
      else ##3 real; I will put as -1 for now
        res1=-1
        res2=-1
        res3=-1
      end

    end
  end
  return res1,res2,res3
end
function cubic3(a::F, b::F, c::F, d::F) where F <: AbstractFloat
  if a == 0.0
    res1,res2 = allrealquadratic(b,c,d)
    res3=Inf
  else
    if d == 0.0
      res1,res2 = allrealquadratic(a,b,c)
      res3=0.0
    else
      A = b/a
      B = c/a
      C = d/a
      Q = (A^2-3B)/9
      R = (2*A^3-9A*B+27C)/54
      p=-3*Q
      q=2*R
      S = -A/3
      if 4p*p*p+27q*q>0  #2Imag 1 real
        m=sqrt((q*q/4)+p*p*p/27)
        res1=cbrt(-q/2+m)+cbrt(-q/2-m)+S
        res2=Inf
        res3=Inf
      else ##3 real; I will put as -1 for now
        res1=-1
        res2=-1
        res3=-1
      end
    end
  end
  return res1,res2,res3
end
function cubic4(a::F, b::F, c::F, d::F) where F <: AbstractFloat#Cardano’s method:
  if a == 0.0
    res1,res2 = allrealquadratic(b,c,d)
    res3=Inf
  else
    if d == 0.0
      res1,res2 = allrealquadratic(a,b,c)
      res3=0.0
    else
          A = b/a
          B = c/a
          C = d/a
          Q = (A^2-3B)/9
          R = (2*A^3-9A*B+27C)/54
          p=-3*Q
          q=2*R
          S = -A/3
          if 4p*p*p+27q*q>0 #2Imag 1 real
                resVect = allrealquadratic(1.0,q,-(p^3)/27)
                z=resVect[1]
              #  @show resVect
                α=cbrt(resVect[1])
                β=-p/(3α)
                res1=α+β+S 
                
               #=  α=cbrt(resVect[2])
                β=-p/(3α)
                res2=α+β+S =#
                res2=Inf
                res3=Inf
          else ##3 real; I will put as -1 for now
            res1=-1
            res2=-1
            res3=-1
          end
      end
  end
  return res1,res2,res3
end
function cubic5(a::F, b::F, c::F, d::F) where F <: AbstractFloat
    _a = one(F) / a
    b, c, d = b * _a, c * _a, d * _a
   #=  m = b < c ? b : c
    m = d < m ? d : m
    m > eps(F) && return typemax(F) # Cauchy bound =#
    _3 = one(F) / 3
    _9 = one(F) / 9
    SQ3 = sqrt(3one(F))
    xₙ = -b * _3
    b²_9 = b * b * _9
    yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)   #eq to 2R
    δ² = muladd(-_3, c, b²_9)                  #eq to Q
    h² = 4δ² * δ² * δ²
    Δ = muladd(yₙ, yₙ, -h²)
    if Δ > 4eps(F) # one real root and two complex roots
      p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
      q = δ² / p
      res1 = xₙ + p + q
      res2=Inf
      res3=Inf
      #z > -eps(F) ? z : typemax(F)
    elseif Δ < -4eps(F) # three real roots
      θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
      δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
      z₁ = 2δ * cos(θ)
      z₂ = muladd(-0.5, z₁, xₙ)
      z₃ = SQ3 * δ * sin(θ)
      res1 = xₙ + z₁
      res2 = z₂ + z₃
      res3 = z₂ - z₃
      #= x = x₁ > -eps(F) ? x₁ : typemax(F)
      x = x₂ > -eps(F) && x₂ < x ? x₂ : x
      x₃ > -eps(F) && x₃ < x ? x₃ : x =#
    else # double repeated real roots
      δ = cbrt(0.5yₙ)
      res1 = xₙ + δ
      res2 = xₙ - 2δ
      res3 =Inf
    #=  x = x₁ > -eps(F) ? x₁ : typemax(F)
      x₂ > -eps(F) && x₂ < x ? x₂ : x =#
    end
    return res1,res2,res3
end
using BenchmarkTools
function closeEnough(sum1::Float64,sum2::Float64) 
    return abs(sum1-sum2)>1e-12
end
#= @show cubic1(-1.33,2.0,1.0,6.0)#1.007151552710968, Inf, Inf)   #  res1= 0.5*(T+U) where T = -sign(R)*cbrt(abs(R)+sqrt(R^2-Q^3))  is wrong
@show cubic2(-1.33,2.0,1.0,6.0) #2.515556238254016, Inf, Inf)
@show cubic3(-1.33,2.0,1.0,6.0) #2.515556238254017, Inf, Inf)
@show cubic4(-1.33,2.0,1.0,6.0) # 2.5155562382540193,Inf, Inf
@show cubic5(-1.33,2.0,1.0,6.0)#2.515556238254016, Inf, Inf) =#

#= @show cubic1(1.33,2.0,1.0,6.0)#(-0.8186606842427812, Inf, Inf)
@show cubic2(1.33,2.0,1.0,6.0) #(-2.1385740543503533, Inf, Inf)
@show cubic3(1.33,2.0,1.0,6.0) #(-2.138574501422444, Inf, Inf)
@show cubic4(1.33,2.0,1.0,6.0) # (-2.1385745013176427, Inf, Inf)
@show cubic5(1.33,2.0,1.0,6.0)#(-2.1385745013176427, Inf, Inf) =#

#= @show cubic1(-1.33,-2.0,1.0,6.0)#(0.9437445952160505, Inf, Inf)
@show cubic2(-1.33,-2.0,1.0,6.0) #(1.3862360576000208, Inf, Inf)
@show cubic3(-1.33,-2.0,1.0,6.0) #(1.386236057600021, Inf, Inf)
@show cubic4(-1.33,-2.0,1.0,6.0) # (1.3862360576000199, Inf, Inf)
@show cubic5(-1.33,-2.0,1.0,6.0)#(1.3862360576000208, Inf, Inf) =#

#= @show cubic1(-1.33,-2.0,-1.0,6.0)#(0.8339202659830023, Inf, Inf)
@show cubic2(-1.33,-2.0,-1.0,6.0) #(1.1665873991339244, Inf, Inf)
@show cubic3(-1.33,-2.0,-1.0,6.0) #(1.1665873987668016, Inf, Inf)
@show cubic4(-1.33,-2.0,-1.0,6.0) # (1.1665861175917736, Inf, Inf)
@show cubic5(-1.33,-2.0,-1.0,6.0)#(1.1665873991339244, Inf, Inf) =#

#= @show cubic1(1.33,12.65,1.0,-6.0)#(-9.379843787886902, 0.6308983740365317, -0.7623327816383512)
@show cubic2(1.33,12.65,1.0,-6.0) #(-1, -1, -1)
@show cubic3(1.33,12.65,1.0,-6.0) #(-1, -1, -1)
@show cubic4(1.33,12.65,1.0,-6.0) # (-1, -1, -1) =#
a=1.33
b=12.65
c=1.0
d=-6.0
  

#@show cubic5(1.33,12.65,1.0,-6.0)# (-9.379843787886902, -0.7623327816383527, 0.6308983740365341)
@btime cubic5(a,b,c,d)# (-9.379843787886902, -0.7623327816383527, 0.6308983740365341)