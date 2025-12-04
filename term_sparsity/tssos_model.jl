module tssos_model
#=
construct the tssos problem model, which is the dual problem of the GMP
=#
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS

src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")

export get_γₜᵗˢˢᵒˢ

###################################################################################
"""
max_{p ∈ R[x]_E }<p,y>
s.t. F-p = σ₀+∑ᵢ₌₁ⁿσᵢxᵢ+τ(x'*x-1)
"""

function modelγₜᵗˢˢᵒˢ(A,n,d,t,F)
    #mosek_setting=mosek_para(1e-8, 1e-8, 1e-8, -1, 1)
    #model = Model(optimizer_with_attributes(Mosek.Optimizer))
    model = Model()
    

    # Define the coefficents of polynomial variable p
    @variable(model, p[1:length(A)])

    # F-p = σ₀+∑ᵢ₌₁ⁿσᵢxᵢ+τ(x'*x-1)-------------------------------------------------------------------
    @polyvar x[1:n]
    h = [x'*x-1]# the equality constraints
    g = [x...] # the inequality constraints

    t₀ = Int(ceil( (d+1)/2 ))
    bracx = vcat([MultivariatePolynomials.monomials(x, i) for i=0:t₀]...) # return all monomials of degree ≤ t₀
    F = bracx'*(F)*bracx
    nonne = F - sum(p.*MultivariatePolynomials.monomials(x, d))

    # QUIET：true then no process information will be output 
    # CS: correlative sparse,"MF" by default (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation)
    # TS: term sparse, "block" by default (maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)
    # SO: sparse order
    #info1 = add_psatz!(model, nonne, x, g, h, t, QUIET=false, CS=true, cliques=[], TS="block", SO=1, Groebnerbasis=false, constrs="con1")
    info1 = add_psatz!(model, nonne, x, g, h, t, QUIET=false, CS=false, cliques=[], TS=false, SO=1, Groebnerbasis=false, constrs="con1")

    # max <p,y>---------------------------------------------------------------------------------------
    subscripts = moments.get_subscripts(d, n)
    Avals = [A[item] for item in subscripts]
    @objective(model, Max, sum(p.*Avals))
    return model, info1
end

function computeγₜᵗˢˢᵒˢ(model)
    set_optimizer(model, Mosek.Optimizer)
    #set_attribute(model, "MSK_IPAR_NUM_THREADS", 1) # if we set the thread of mosek to be 1, it seems that the tssos model will not become faster in the julia mode
    
    optimize!(model)

    println("sol_time: ", solve_time(model))
    println("status: ",termination_status(model))
    println("Primal: ", primal_status(model)) # primal_status(model) >>FEASIBLE_POINT::ResultStatusCode = 1 indicates that it find a feasible point
    println("Dual: ", dual_status(model))
    println("Ojbective: ", objective_value(model))
    return model
end

"""Computes one of the many flavour of γₜᵗˢˢᵒˢ"""
function get_γₜᵗˢˢᵒˢ(A,n,d,t,F)
    mod, info1 = modelγₜᵗˢˢᵒˢ(A,n,d,t,F)
   
    γₜᵗˢˢᵒˢ = computeγₜᵗˢˢᵒˢ(mod)
    ex_mom = -dual(constraint_by_name(mod, "con1")) # the value of tms Y
    return γₜᵗˢˢᵒˢ, ex_mom, info1
end

end