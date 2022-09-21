module qssgenerated
using RuntimeGeneratedFunctions
using StaticArrays
using SymEngine
using ExprTools
using MacroTools: postwalk, prewalk, @capture
import Base.:-
import Base.:+
import Base.:*
using Plots: plot!


RuntimeGeneratedFunctions.init(@__MODULE__)


#this section belongs to taylorseries subcomponent
import Base: ==, +, -, *, /, ^

import Base: iterate, size, eachindex, firstindex, lastindex,
    eltype, length, getindex, setindex!, axes, copyto!

import Base:  sqrt, exp, log, sin, cos, sincos, tan,
    asin, acos, atan, sinh, cosh, tanh, atanh, asinh, acosh,
    zero, one, zeros, ones, isinf, isnan, iszero,
    convert, promote_rule, promote, show,
    real, imag, conj, adjoint,
    rem, mod, mod2pi, abs, abs2,
   
    power_by_squaring,
    rtoldefault, isfinite, isapprox, rad2deg, deg2rad
    #end of taylorseries section of imports


    # list of public (API) to the user, not between files as those are linked as if in one file
    export SimSettings,QSS_Problem,QSS_Solve ,  qss1,qss2,qss3,liqss1,liqss2,saveat,plotSol,evaluateSol

    export  @NLodeProblem,save_prob_to_model,QSS_Solve_from_model

    export Taylor0,mulT,mulTT,createT,addsub,negateT,subsub,subadd,subT,addT,muladdT,mulsub,divT


    #include section of ts subcomponent
    include("ownTaylor/parameters.jl")  
    include("ownTaylor/constructors.jl") 
    include("ownTaylor/conversion.jl")        
    include("ownTaylor/arithmetic.jl")
    include("ownTaylor/arithmeticT.jl")
    include("ownTaylor/functions.jl")
    include("ownTaylor/functionsT.jl")
    include("ownTaylor/power.jl")
    include("ownTaylor/powerT.jl")
    include("ownTaylor/auxiliary.jl") 
    include("ownTaylor/calculus.jl")    
    include("ownTaylor/other_functions.jl")
    include("ownTaylor/evaluate.jl")
       
       

    #Utils
    include("Utils/SimUtils.jl") 
    #QSSFamily/common/
    include("Common/Solution.jl")
    include("Common/QSSNLProblemHelper.jl")
    include("Common/QSSNLProblem.jl")

    include("Common/QSS_data.jl")
    include("Common/Scheduler.jl")
    include("Common/QSS_modelworkshop.jl")


    include("NL_integrators/NL_QSS_Integrator.jl")
    include("NL_integrators/NL_LiQSS_Integrator.jl")
    include("Quantizers/QSS_quantizer.jl")
    include("Quantizers/LiQSS_quantizer.jl")


end # module
