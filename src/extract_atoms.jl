module extract_atoms

using TypedPolynomials
using MultivariateMoments#, MultivariatePolynomials
using HomotopyContinuation

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")

export get_atoms

"""takes the moment df and extracts: n, t, and moment values"""
proc_mom(df::Matrix{Any}) = typeof(df[1]) == Vector{Any} ? proc_sparse_mom(df) : proc_dense_mom(df)
function proc_sparse_mom(df::Matrix{Any})
    m_i, m_v = get_inds_vals(df)
    nc = maximum([m[1] for m in m_i]) # the number of cliques
    n_s = [length([m[2] for m in m_i if m[1] == c][1]) for c in 1:nc] # the number of variables of each clique
    t = div(maximum([m[2][1] for m in m_i if m[1]==1]), 2)
    sz = length(m_v)
    val_s = [[m_v[k] for k in 1:sz if m_i[k][1] == c] for c in 1:nc]
    return n_s, t, val_s
end
function proc_dense_mom(df::Matrix{Any})
    m_i, m_v = get_inds_vals(df) # Obtain the index and corresponding value of tms
    n = length(m_i[1]) # obtain the dimension of tensor
    t = div(maximum(maximum([m for m in m_i])),2) # obtain the order 2*t of tms
    # Note here that val is need to be rewritten as below,
    # otherwhise, the type of m_v itself is Vector{Any}
    # which will lead error using LinearAlgebra.svd(getmat(Moment)), when extractatoms
    val = [item for item in m_v]
    return n, t, val
end
get_inds_vals(df) = df[:,1],df[:,2]

"""Builds moment matrix and extracts atoms"""
get_atoms(n::Vector{Int64}, t, mom, tol=1e-4, hom=true) = [get_atoms(n[j], t, mom[j], tol, hom) for j in 1:length(n)]
function get_atoms(n::Int64, t, mom, tol=1e-4, hom=true)
    M = get_mom_mat(n,t,mom)
    if hom
        solver = HomotopyContinuation.SemialgebraicSetsHCSolver(; compile = true)
        return extractatoms(M,tol,solver) # Call the atomic_measure in extract.jl of MultivarMultivariateMoment
    else
        return extractatoms(M,tol) # the default solver of SemialgebraicSets is used which currently computes the Gröbner basis, 
                                   # then the multiplication matrices and then the Schur decomposition of a random combination of these matrices. 
                                   # For floating point arithmetics, homotopy continuation is recommended as it is more numerically stable than Gröbner basis computation.
    end
end

function get_mom_mat(n::Int64, t::Int64, mom)
    x, μ = get_psudomeasure(n,t,mom)
    B = map(a -> prod(x.^a), moments.make_mon_expo(n,t))
    # return a moment matrix (μ(BBᵀ))
    return MultivariateMoments.moment_matrix(μ, B) # Note the order in this moment matrix is from x[n] to x[1]
end

function get_psudomeasure(n,t,mom_vals)
    x, monomials = get_psudome_set_up(n,t) # the n dimensional variable x of a polynomial, all monomials of degree ≤ 2*t
    return x, MultivariateMoments.measure(mom_vals[1:length(monomials)],monomials)
end

function get_psudome_set_up(n,t)
    @HomotopyContinuation.polyvar x[1:n]
    monomials_expos = moments.make_mon_expo(n,2*t)
    monomials = map(a -> prod(x.^a), monomials_expos)
    return x, monomials
end

"""Extracts the centers and weights of the atom objects"""
ext_centers_weights(extract) = [a.center for a in extract.atoms], [a.weight for a in extract.atoms]
function ext_centers_weights(ext,A,n,d)
    mc = moments.get_maximal_cliques(A,n,d)
    p = length(mc)
    cents = [[[i in mc[j] ? popfirst!(a.center) : 0 for i in 1:n] for a in ext[j].atoms] for j in 1:p]
    weights = [[a.weight for a in ext[j].atoms] for j in 1:p]
    return cents, weights
end

"""(atoms,weights) -> recovered tms corresponding to cp_tensor"""
function recon_ten_id(c,w,d)
    n = size(c[1],1)
    A = Dict()
    subscripts = moments.get_subscripts(d,n)
    for item in subscripts
        A[item] = sum([ w[j]*prod([c[j][i] for i in item]) for j in 1:size(c,1) ])
    end
    return A
end
function recon_ten_sp(c,w,d)
    n = size(c[1][1],1)
    A = Dict()
    subscripts = moments.get_subscripts(d,n)
    for item in subscripts
        A[item] = sum([ sum([ w[k][j] * prod([c[k][j][i] for i in item]) for j in 1:size(w[k],1) ]) for k in 1:size(w,1)])
    end
    return A
end


end