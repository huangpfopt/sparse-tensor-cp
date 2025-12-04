module moments
proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"tensorGraphs.jl")
using .tensorGraphs

export eᵢ,
       make_mon_expo,
       get_maximal_cliques,
       get_monomial_cliques



"""The standard basis vector  eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]

# make_mon_expo----------------------------------------------------------------
"""[x]≤ₜ := [xᵅ for all α ∈ Nⁿ≤ₜ] or [x]₌ₜ:= [xᵅ for all α ∈ Nⁿ₌ₜ]"""
#the generated power is from low to high, for the same power, the index is given by lexicographical order（i.e. Sort first by the largest first component, then by the largest second component.）
function make_mon_expo(n::Int,t::Int; isle::Bool = true) # behind ; is the key parameter, when it is not given, use defaul value, to pass a value need to specify isle=xxx
    t == 0 ? (return[eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i in 1:n, m in tmp],:,1)
    return unique(isle ? vcat(tmp, M_vec) : M_vec)
    # unique is necessary, for example, [1,0],[0,1] add one both can obtain [1,1]
end
"""[x]≤ₜ[x]≤ₜᵀ or [x]₌ₜ[x]₌ₜᵀ"""
make_mon_expo(mom₁::Vector{Vector{Int64}},mom₂::Vector{Vector{Int64}}) = [a+b for a in mom₁, b in mom₂]
make_mon_expo(mom::Vector{Vector{Int64}}) = make_mon_expo(mom,mom)

# get_monomial_cliques---------------------------------------------------------
"""[x(Vₖ)]ₜ"""
function get_monomial_cliques(t::Int,mc::Vector{Vector{Int64}},n::Int)
    Iks = [get_mon_clique(n,t,c) for c in mc]
    Iks_comp = setdiff(make_mon_expo(n,t), union(Iks...))
    return [Iks..., Iks_comp]
end
"""[x(Vₖ)]ₜ[x(Vₖ)]ₜᵀ"""
get_monomial_cliques(t::Tuple{Int,Int},mc,n) = [make_mon_expo(m,m) for m in get_monomial_cliques(t[1],mc,n)]
get_mon_clique(n,t,c) = map(v->embed(v,c,n), make_mon_expo(size(c,1),t))
# v is the index of a tms in a clique, c is the set of nodes for some clique
embed(v,c,n) = [i in c ? popfirst!(v) : 0 for i in 1:n]

## cliques --------------------------------------------------------------------
function get_maximal_cliques(A, n, d, isps=false)
    if isps
        # ps case not implemented yet
    else
        return tensorGraphs.maximal_cliques(A,n,d)
    end
end

## Utilities ------------------------------------------------------------------
"""get all nonrepeat subscripts of of d-order n-dimensional symmetric tensor"""
function get_subscripts(d::Int64, n::Int64,is_cp=true)
    if is_cp
        len = binomial(n+d-1,d)
        subscripts = Vector{Int64}[] #the type of each element in the list is Vector{Int64}
        v = [1 for _ in 1:d]
        for _ in 1:len
            push!(subscripts,v)
            v = v[:] # when v change below, it will not influence those v, which have been pushed into subscripts
            for i in d:-1:1
                if v[i] < n
                    v[i] += 1
                    for j in i+1:d
                        v[j] = v[i]
                    end
                    break
                end
                # only when v[i] == n, the index at the i-1th position will be added by 1
            end
        end
        return subscripts
    else
        return # ps case not implement yet
    end
end

"""get all [i,j,...,l] s.t. Aᵢⱼ⋯ₗ ≠ 0"""
function get_nonzero_entries(A)
    return [item for item in keys(A) if A[item] != 0]
end

end