module cp_model
#=
construct the cp problem model
=#
using JuMP
using MosekTools

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
using .moments; const mom = moments

export get_ξₜᶜᵖ

#########################################################################################

"""
    min L(F) where F is a random SOS function
    s.t L(x∘⋯∘x) = A, i.e. L(xᵢxⱼ⋯xₗ)=Aᵢⱼ⋯ₗ
        L([x]ₜ[x]ₜᵀ) ⪰ 0
        L( (xᵀx-1)[x]ₜ₋₁[x]ₜ₋₁ᵀ ) = 0
        L( xᵢ[x]ₜ₋₁[x]ₜ₋₁ᵀ ) ⪰ 0, i ∈ [n]
"""

# Construct the model without ideal sparsity, A is the corresponding tms of given tensor, t is the relaxation order，F is the coefficent matrix of random generated SOS 
function modelξₜᶜᵖⁱᵈ(A,n,d,t,F)
    #n is dimension of tensor, d is order of tensor
    model = Model()

    # Define tms variable Y = L([x]₂ₜ)----------------------------------------------
    # make_mon_expo(n,2*t) return the list of all power index α with |α|≤2t
    # the index of Y[~] is the genereted α
    @variable(model, Y[mom.make_mon_expo(n,2*t)])

    # L([x]ₜ[x]ₜᵀ) ⪰ 0-----------------------------------------------------
    # make_mom_mat_con(A,t,Y) will construct the moment matrix Mₜ(Y)
    @constraint(model, make_mom_mat_con(n,t,Y) in PSDCone())

    # L(x∘⋯∘x) = A--------------------------------------------------------
    for (c,v) in make_d_ord_mom_con(A,n,Y) @constraint(model, c == v) end #the bracket of (c,v) is neccesary，tuple will be given to c,v automatically

    # L( (xᵀx-1)[x]ₜ₋₁[x]ₜ₋₁ᵀ ) = 0 is eqaual to L( (xᵀx-1)[x]₂ₜ₋₂ )=0-------------
    @constraint(model,make_norm_con(n,t,Y) .== 0)

    # L( xᵢ[x]ₜ₋₁[x]ₜ₋₁ᵀ ) ⪰ 0, i ∈ [n]-------------------------------------
    for c in make_loc_con(n,t,Y)
        t == 1 ? @constraint(model, c[1] ≥ 0) : @constraint(model, c in PSDCone())
    end

    # min L(F)--------------------------------------------------------------
    @objective(model, Min, make_obj(n,d,Y,F))
    return model
end

"""Construnct the moment matrix Mₜ(Y)"""
# a_to_yₐ(Y,a) will map the indexes in a to corresponding Y[~]
# mom.make_mom_expo(mom::Vector{Vector{Int64}}) return the moment matrix
make_mom_mat_con(n,t,Y) = a_to_yₐ(Y, mom.make_mon_expo(mom.make_mon_expo(n,t)) )

"""Obtain the one to one map between L(x∘⋯∘x) and the elements of A"""
function make_d_ord_mom_con(A,n,Y)
    return [(Y[sub2pow(ij,n)], A[ij]) for ij in keys(A)]
end

"""Construct L( (xᵀx-1)[x]ₜ₋₁[x]ₜ₋₁ᵀ ), L( (xᵀx-1)[x]₂ₜ₋₂ )"""
function make_norm_con(n,t,Y)
    # Construct the index of the tms with degree 2t-2
    mon₂ₜ₋₂ = mom.make_mon_expo(n,2t-2)
    Lₓₜₓ₋₁(α) = sum([Y[α+2*mom.eᵢ(n,i)] for i in 1:n])-Y[α]
    return map(α -> Lₓₜₓ₋₁(α), mon₂ₜ₋₂)
end

"""Construct L( xᵢ[x]ₜ₋₁[x]ₜ₋₁ᵀ )"""
function make_loc_con(n,t,Y)
    momₜ₋₁ = mom.make_mon_expo(mom.make_mon_expo(n,t-1))
    Lₓᵢ(i) = map(α -> Y[α+mom.eᵢ(n,i)],momₜ₋₁)
    return [Lₓᵢ(i) for i in 1:n]
end

"""<F,Y>"""
function make_obj(n,d,Y,F)
    t₀ = Int(ceil( (d+1)/2 ))
    momₜ₀ = mom.make_mon_expo(mom.make_mon_expo(n,t₀))
    return sum(F.*a_to_yₐ(Y, momₜ₀))
end

####################################################################################

"""
    min ∑ₖ Lᵥₖ(F) where F is a random SOS function
    s.t ∑ₖ Lᵥₖ(x∘⋯∘x) = A, 
        Lᵥₖ([x]ₜ[x]ₜᵀ) ⪰ 0, k ∈ [p]
        Lᵥₖ( (xᵀx-1)|_ᵥₖ[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ) = 0, k ∈ [p]
        Lᵥₖ( xᵢ[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ) ⪰ 0, i ∈ Vₖ, k ∈ [p]
"""

function modelξₜᶜᵖˢᵖ(A,n,d,t,F)
    model = Model()

    # Obtain the maximal cliques of the tensor
    # that is the index set of variables that may be nonzero in the same decomposition vector
    mc = get_maximal_cliques(A,n,d)
    # generate the index α |α|≤2t for tms in each clique
    spar_inds = make_spar_inds(mc,t)
    # Define tms variable Y，L([x(Vₖ)]₂ₜ), k ∈ [p]---------------------------------------------
    # spar_inds is the index of Y 
    @variable(model, Y[spar_inds])

    # Lᵥₖ([x]ₜ[x]ₜᵀ) ⪰ 0, k ∈ [p]-----------------------------------------------------
    # make_spar_mom_mat_con(A,t,Y) will construct the moment matrix M(Yₖ), k ∈ [p] corresponding to each cliques
    for m in make_spar_mom_mat_con(A,n,d,t,Y) @constraint(model, m in PSDCone()) end

    # ∑ₖ Lᵥₖ(x∘⋯∘x) = A---------------------------------------------------------------
    for (c,v) in make_spar_d_ord_mom_con(A,n,d,Y) @constraint(model, c == v) end

    # Lᵥₖ( (xᵀx-1)|_ᵥₖ[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ) = 0, k ∈ [p]-----------------------------------
    for c in make_spar_norm_con(A,n,d,t,Y) @constraint(model, c .== 0) end

    # Lᵥₖ( xᵢ[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ) ⪰ 0, i ∈ Vₖ, k ∈ [p]------------------------------------
    for c in make_spar_loc_con(A,n,d,t,Y)
        t == 1 ? @constraint(model, c[1] ≥ 0) : @constraint(model, c in PSDCone())
    end

    # min ∑ₖ Lᵥₖ(F)---------------------------------------------------------------------
    @objective(model, Min, make_spar_obj(A,n,d,Y,F))
    return model
end

make_spar_inds(mc,t) = [ [k,m] for k in 1:size(mc,1) for m in mom.make_mon_expo(size(mc[k],1), 2*t)]

"""Lᵥₖ([x]ₜ[x]ₜᵀ) ⪰ 0, k ∈ [p]"""
function make_spar_mom_mat_con(A,n,d,t,Y)
    mc = get_maximal_cliques(A,n,d)
    p = size(mc,1)
    # Obtain the index of Mₜ(Yₖ) corresponding to each clique, and embed these index to the index of full n-dimensional vectors
    Iᵏs = mom.get_monomial_cliques((t,t),mc,n)[1:p]
    Iᵏs_cut = [map(a-> [k, a[mc[k]] ], Iᵏs[k]) for k in 1:p]
    return [a_to_yₐ(Y, Iᵏ) for Iᵏ in Iᵏs_cut]
end

"""∑ₖ Lᵥₖ(x∘⋯∘x) = A"""
function make_spar_d_ord_mom_con(A,n,d,Y)
    mc = get_maximal_cliques(A,n,d)
    p = size(mc,1)
    nze = mom.get_nonzero_entries(A)  
    x_d = map(ij -> sub2pow(ij, n), nze)
    K = sum([map(α -> isinmc(α,mc[k]) ? Y[[k,α[mc[k]] ]] : 0, x_d) for k in 1:p])
    A_value = [A[ij] for ij in nze]
    return [(K[i], A_value[i]) for i in 1:size(nze,1)]
end
isinmc(α,mck) = all([j in mck for j in findall(α.>0)])

"""Lᵥₖ( (xᵀx-1)|_ᵥₖ[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ) = 0, k ∈ [p]"""
function make_spar_norm_con(A,n,d,t,Y)
    mc = get_maximal_cliques(A,n,d)
    p = size(mc,1)
    #Construct the index of tms of degree 2t-2 corresponding to each clique, and embed these index to the index of full n-dimensional vectors
    Iᵏs = mom.get_monomial_cliques(2t-2,mc,n)
    Lᵏₓₜₓ₋₁(α,k) = sum([Y[[k,(α+2*mom.eᵢ(n,i))[mc[k]] ]] for i in mc[k]]) - Y[[k,α[mc[k]] ]]
    return [[Lᵏₓₜₓ₋₁(α,k) for α in Iᵏs[k]] for k in 1:p]
end

"""Lᵥₖ( xᵢ[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ) ⪰ 0, i ∈ Vₖ, k ∈ [p]"""
function make_spar_loc_con(A,n,d,t,Y)
    mc = get_maximal_cliques(A,n,d)
    p = size(mc,1)
    #Construct the index of tms of degree t-1 corresponding to each clique, and embed these index to the index of full n-dimensional vectors
    Iᵏs = mom.get_monomial_cliques(t-1,mc,n)
    Lᵏₓᵢ(i,α,k) = Y[[k,(α+mom.eᵢ(n,i))[mc[k]] ]]
    arr = [[[Lᵏₓᵢ(i,α+β,k) for α in Iᵏs[k], β in Iᵏs[k]]
                           for i in mc[k]]
                           for k in 1:p]
    return cat(arr...,dims=1)
end

"""min ∑ₖ Lᵥₖ(F)"""
function make_spar_obj(A,n,d,Y,F)
    mc = get_maximal_cliques(A,n,d)
    p = size(mc,1)
    t₀ = Int(ceil( (d+1)/2 ))
    # the index of moment matrix Mₜ₀(Y)
    momₜ₀ = mom.make_mon_expo(mom.make_mon_expo(n,t₀))
    # the index of the t₀-th moment matrix Mₜ₀(Yₖ) corresponding to each clique, and embed these index to the index of full n-dimensional vectors
    Iᵏsₜ₀ = mom.get_monomial_cliques((t₀,t₀),mc,n)[1:p]
    K = sum([map(α -> (α in Iᵏsₜ₀[k]) ? Y[[k,α[mc[k]] ]] : 0, momₜ₀) for k in 1:p])
    return sum(F.*K)
end

####################################################################################
"""This is where the ξₜᶜᵖ is calculated for tms A"""
function computeξₜᶜᵖ(model)
    set_optimizer(model, Mosek.Optimizer) # Use Mosek as the sdp solver
    #set_attribute(model, "MSK_IPAR_NUM_THREADS", 1) #Set the thread of mosek to be 1, this setting has no impact on the result, guess might be that the sparse case do not need multi-thread
    optimize!(model)
    #https://jump.dev/JuMP.jl/stable/manual/solutions/ for JuMP solution information
    println("sol_time: ", solve_time(model)) #returns the solve time in wall-clock seconds reported by the solver
    println("status: ",termination_status(model)) #  is_solved_and_feasible(model)
    println("Primal: ", primal_status(model)) # primal_status(model) >>FEASIBLE_POINT::ResultStatusCode = 1 indicates that it find a feasible point
    println("Dual: ", dual_status(model))
    println("Ojbective: ", objective_value(model))
    return model
end

"""Computes one of the many flavour of ξₜᶜᵖ"""
function get_ξₜᶜᵖ(A,n,d,t,F,flavour)
    if contains(flavour,"id")
        mod = modelξₜᶜᵖⁱᵈ(A,n,d,t,F)
    elseif contains(flavour,"sp")
        mod = modelξₜᶜᵖˢᵖ(A,n,d,t,F)
    else
        error("Incorrect model specification")
    end
    ξₜᶜᵖ = computeξₜᶜᵖ(mod)
    ex_mom = extract_moments(ξₜᶜᵖ) # the value of tms Y
    return ξₜᶜᵖ, ex_mom
end

"""Returns L(xᵅ) for α ∈ Nⁿ₂ₜ"""
function extract_moments(model)
    Y = model[:Y] # access the regestered variable Y of model
    indices = model[:Y].axes[1] # {a ∈ model ⊆ Nⁿ₂ₜ or Nᵛ₂ₜ}
    values = JuMP.value.(Y).data # {yₐ ∈ model}
    return hcat(indices, values)
end

######################### Utilities #################################################
"""(L,[α]ᵢⱼ) -> [L(xᵅ)]ᵢⱼ"""
a_to_yₐ(Y, i_arr) = map(a -> Y[a], i_arr)

"""subscript of tensor [i,j,...,l] -> power α"""
function sub2pow(subscript,n)
    alpha = [0 for _ in 1:n]
    for i in subscript
        alpha[i] += 1
    end
    return alpha
end

end