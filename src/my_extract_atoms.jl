module my_extract_atoms

using LinearAlgebra
using RowEchelon

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
get_atoms(n::Vector{Int64}, t::Vector{Int64}, mom, tol=1e-4) = [get_atoms(n[j], t[j], mom[j], tol) for j in 1:length(n)]
# In our current examples, the result of using t::Vector{Int64} and t seems almost the same
get_atoms(n::Vector{Int64}, t, mom, tol=1e-4) = [get_atoms(n[j], t, mom[j], tol) for j in 1:length(n)]
function get_atoms(n::Int64, t, mom, tol=1e-4)
    # Obtain the t-th order moment matrix
    M = get_mom_mat(n,t,mom)
    
    #Step1: COmpute the Choleskey decomposition of M=UUᵀ 
    F = svd(M) # M = F.U*Diagonal(F.S)*F.Vt, the singular values are sorted in decreasing order
    #F = eigen(M,sortby=-) # the result seems is similar to using svd M = F.vectors*Diagonal(F.values)*F.vectors', "sortby=-" to sort eigenvalues in decreasing order
    #nM = F.S[1] # norm of M 
    
    r = something(findfirst(σ -> σ < tol, F.S), 0) # If the return value of findfirst is not empty, return the value, otherwise something will return 0
    #r = something(findfirst(σ -> σ < tol, F.values), 0)
    if iszero(r)
        #cM = tol
        r = length(F.S)
        #r = length(F.values)
    else
        #cM = F.S[r] / nM
        r -= 1
    end
    
    #r = rank(M)
    U = F.U[:,1:r] * Diagonal(sqrt.(F.S[1:r]))
    #U = F.vectors[:,1:r] * Diagonal(sqrt.(F.values[1:r]))

    # Step2: reduce U to column echelon form
    # Row-reduced echelon form of the transpose of U
    U, pivots = RowEchelon.rref_with_pivots!(Matrix(U'),tol)
    U = Matrix(U')

    # Step3: extract multiplication matrix for each variable
    Iₜ = moments.make_mon_expo(n, t)
    BasisPow = Iₜ[pivots]
    Ns = []
    for i=1:n
        rows = [findfirst(a -> a == (BasisPow[j]+moments.eᵢ(n,i)), Iₜ) for j in 1:r]
        push!(Ns, U[rows,:])
    end

    # random combination of multiplication matrices
    coef = rand(n)
    coef = coef/sum(coef)
    combNs = sum([coef[i]*Ns[i] for i in 1:n])

    # Step4: perform ordered schur decomposition of combNs, eigenvalues sorted in increasing order along the diagonal of T
    #F = schur(combNs)
    #position = sortperm(diag(F.T)) 
    #Q = F.Z[:,position] # wrong way
    (Q,T) = orderschur(combNs)

    # Step5: retrive optimal vectors 
    cents = [[Q[:,j]'*Ns[i]*Q[:,j] for i in 1:n] for j in 1:r]

    # compute the weights
    I₂ₜ = moments.make_mon_expo(n, 2*t)
    A = zeros(binomial(n+2t,2t), r)
    for i in 1:r 
        A[:,i] = map(a->prod(cents[i].^a),I₂ₜ)
    end
    weights = A \ (mom[1:length(I₂ₜ)])
    return [cents, weights]
end

get_atoms_robust(n::Vector{Int64}, t, mom, tol=1e-4) = [get_atoms_robust(n[j], t, mom[j], tol) for j in 1:length(n)]
function get_atoms_robust(n::Int64, t, mom, tol=1e-4) # this function was written based on the old version (v1.3.4?) of Tssos.extract_solutions_robust, but later Tssos was updated on GitHub.
    # Obtain the t-th order moment matrix
    M = get_mom_mat(n,t,mom)
    Iₜ = moments.make_mon_expo(n, t)
    ls = binomial(n+t-1, n)
    N = Vector{Matrix{Float64}}(undef, n)
    for i = 1:n
        N[i] = zeros(Float64, ls, ls)
        temp = moments.eᵢ(n,i)
        for j = 1:ls, k = j:ls
            #loc = bfind_to(basis, size(basis, 2), basis[:,k] + temp, n)
            loc = findfirst(a -> a == (Iₜ[k]+moments.eᵢ(n,i)), Iₜ) 
            N[i][j,k] = M[j, loc]
        end
        N[i] = Symmetric(N[i],:U) # Construct a symmetric matrix via the up trianglar part of N[i]
    end
    F = svd(M[1:ls, 1:ls])
    S = F.S[F.S .>= tol]
    S = sqrt.(S).^(-1)
    for i = 1:n
        N[i] = Diagonal(S)*F.Vt[1:length(S),:]*N[i]*F.U[:,1:length(S)]*Diagonal(S)
    end
    rands = rand(n)
    rands = rands/sum(rands)
    L = schur(sum(rands[i]*N[i] for i in 1:n)).Z
    cents = [[L[:,i]'*N[j]*L[:,i] for j = 1:n] for i = 1:length(S)]

    # compute the weights
    I₂ₜ = moments.make_mon_expo(n, 2*t)
    A = zeros(binomial(n+2t,2t), length(S))
    for i in 1:length(S) 
        A[:,i] = map(a->prod(cents[i].^a),I₂ₜ)
    end
    weights = A \ mom
    return [cents, weights]
end

function rank_check(n::Vector{Int64}, t, d, mom, tol=1e-6) 
    rs = [rank_check(n[j], t, d, mom[j], tol) for j in 1:length(n)]
    check = all([item[1] for item in rs])
    return check,rs
end
function rank_check(n::Int64, t, d, mom, tol=1e-6)
    M = get_mom_mat(n,t,mom)
    rs = Int64[]
    oldrank = 0
    check = false
    #for j = 1:t
    for j = Int(ceil(d/2))-1:t 
        sz = binomial(n+j,j)
        
        F = svd(M[1:sz,1:sz]) # LinearAlgebra.svd
        r = something(findfirst(σ -> σ < tol, F.S), 0) # If the return value of findfirst is not empty, return the value, otherwise something will return 0
        if iszero(r)
            r = length(F.S)
        else
            r -= 1
        end
        
        #r = rank(M[1:sz,1:sz])
        push!(rs, r)
        if r == oldrank
            check = true
            break
        end
        oldrank = r
    end
    return check, rs
end

function get_mom_mat(n::Int64, t::Int64, mom)
    μ = get_psudomeasure(n,t,mom)
    # the index corresponding to the moment matrix Mₜ
    Iₜ = moments.make_mon_expo(moments.make_mon_expo(n,t))
    # return a moment matrix (μ(BBᵀ))
    return map(α->μ[α], Iₜ) 
end

function get_psudomeasure(n,t,mom_vals)
    # an n-dimensional vector of the indexes for all degree ≤ 2*t
    I₂ₜ = moments.make_mon_expo(n, 2*t)
    sz = binomial(n+2t,2t)
    μ = Dict()
    for i in 1:sz 
        μ[I₂ₜ[i]] = mom_vals[i]
    end
    return μ
end

"""Extracts the centers and weights of the atom objects"""
ext_centers_weights(extract) = extract[1], extract[2]
function ext_centers_weights(ext,A,n,d)
    mc = moments.get_maximal_cliques(A,n,d)
    p = length(mc)
    cents = [[[i in mc[j] ? popfirst!(a) : 0 for i in 1:n] for a in ext[j][1]] for j in 1:p]
    weights = [ext[j][2] for j in 1:p]
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


function orderschur(combNs)
    """
    Algorithm 7.6.1 in [Golub. Matrix computation]
    """
    (T,Q) = schur(combNs)
    n = size(combNs,1)

    swap = true
    while swap #bubble sort
        swap = false
        for k=1:n-1
            if T[k,k] > T[k+1,k+1]
                swap = true # the swap happend
                # Givens rotation
                (c,s) =givens(T[k,k+1],T[k+1,k+1]-T[k,k])
                T[k:k+1,k:n] = [c s;-s c]'*T[k:k+1,k:n]
                T[1:k+1,k:k+1] = T[1:k+1,k:k+1]*[c s; -s c]
                Q[1:n,k:k+1] = Q[1:n,k:k+1]*[c s; -s c]
            end
        end
    end
    return Q,T
end

function givens(a,b)
    """
    Algorithm 5.1.3 in [Golub. Matrix computation]
    """
    if b==0
        c=1
        s=0
    else
        if abs(b) > abs(a)
            t=-a/b; s=1/sqrt(1+t^2); c=s*t 
        else
            t = -b/a; c=1/sqrt(1+t^2); s=c*t;
        end
    end
    return c,s
end

end