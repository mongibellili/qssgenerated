using LinearAlgebra
#=
compare: case of 3 reals
cubics based on the depressed cubic  x^3+px+q=0 ;or x^3- 3Qx + 2R = 0
cubic1: if R^2 < Q^3  res1 = -2*sqrt(Q)*cos(ϕ/3)+S where  ϕ = acos(R/sqrt(Q^3)) ;res2 ϕ+2π res2  ϕ-2π
cubic5: if R^2 < Q^3  res1=2√abs(Q) * cos(θ)  where θ = abs(2R) < eps(F) ? 0.5π/3 : atan(√abs(Δ) / abs(q))/3  ;Δ=4(R^2-Q^3))
=#
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
                    h3=Q^3
                   #=  @show 2R
                    @show abs(2R) < eps(1.0) =#
                    S = -A/3
                    if R^2 < Q^3   #3 real
                            P = -2*sqrt(Q)
                            ϕ = acos(R/sqrt(h3))
                           # res1 = P*cos(ϕ/3)+S   
                            res1 = muladd(P,cos(ϕ/3),S )  
                                            
                            res2 = muladd(P,cos((ϕ+2π)/3),S )     
                           # res2 = P*cos((ϕ+2π)/3)+S                       
                            res3 =muladd(P,cos((ϕ-2π)/3),S)
                           # res3 = P*cos((ϕ-2π)/3)+S
                            
                    else   #2Imag 1 real
                            
                            res1=Inf
                            res2=Inf
                            res3=Inf                                           
                    end
                end
            end
  return res1,res2,res3
end #1real2Imag wrong for this algorithm

function cubic5(a::F, b::F, c::F, d::F) where F <: AbstractFloat
    #yₙ =2R = q for other algs
    #δ² = -Q = p/3
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
    if Δ < -4eps(F) # three real roots
      θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3  # 62.157 ns                             # acos(-yₙ / √h²)  
      #θ = acos(yₙ/sqrt(h²))/3   #67.703 ns 
      δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
      z₁ = 2δ * cos(θ)
      z₂ = muladd(-0.5, z₁, xₙ)
      z₃ = sqrt(3) * δ * sin(θ)
      res1 = xₙ + z₁ 
      res2 = z₂ + z₃
      res3 = z₂ - z₃ 
   #=    res1 =  2δ * cos(θ) + xₙ
      res2 =  muladd(-0.5, 2δ * cos(θ), xₙ) + sqrt(3) * δ * sin(θ)
      res3 =  muladd(-0.5, 2δ * cos(θ), xₙ) - sqrt(3) * δ * sin(θ) =#
      #= x = x₁ > -eps(F) ? x₁ : typemax(F)
      x = x₂ > -eps(F) && x₂ < x ? x₂ : x
      x₃ > -eps(F) && x₃ < x ? x₃ : x =#
    else
      res1=-1.0
      res2=-1.0
      res3=-1.0
    #=  x = x₁ > -eps(F) ? x₁ : typemax(F)
      x₂ > -eps(F) && x₂ < x ? x₂ : x =#
    end
    return res1,res2,res3
end
using BenchmarkTools
function closeEnough(sum1::Float64,sum2::Float64) 
    return abs(sum1-sum2)>1e-12
end
#= @show cubic1(1.33,12.65,1.0,-6.0)#(-9.379843787886902, 0.6308983740365317, -0.7623327816383512)
@show cubic5(1.33,12.65,1.0,-6.0)# (-9.379843787886902, -0.7623327816383527, 0.6308983740365341) =#

#= @btime cubic1(1.33,12.65,1.0,-6.0)#(-9.379843787886902, 0.6308983740365317, -0.7623327816383512)
@btime cubic5(1.33,12.65,1.0,-6.0)# (-9.379843787886902, -0.7623327816383527, 0.6308983740365341) =#


#@btime allrealquadratic(1.33,1.0,-6.0)#(-2.5329305695911355, 1.7810508703430152)#1.647 ns (0 allocations: 0 bytes)