module cp_tensors
using Random
Random.seed!(1234) #generate the random sequence that can be repeated

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"matrix_IO.jl")
include(proj_dir*"moments.jl")

export generate_lit_cp_tens,
       generate_lit_non_cp_tens,
       generate_random_cp_tens

## Examples
### Literature Examples
function generate_lit_cp_tens()
    lit_cp_tens = Dict()
    #17Fan, P7
    indices = [[1,1,1],[1,1,2],[1,1,3],[1,2,2],[1,2,3],[1,3,3],
    [2,2,2],[2,2,3],[2,3,3],[3,3,3]]
    values = [2,1,1,1,0,1,2,0,0,1]
    lit_cp_tens["ex01_n3_d3_17Fan-p7"] = hcat(indices,values)

    #17Fan, Example 4.3
    indices = [[2,2,2],[2,2,3],[2,2,4],[2,2,5],[2,2,8],[2,3,3],
    [2,3,8],[2,4,4],[2,4,5],[2,5,5],[2,8,8],[3,3,3],
    [3,3,4],[3,3,5],[3,3,7],[3,3,8],[3,4,4],[3,4,5],
    [3,5,5],[3,7,7],[3,7,8],[3,8,8],[4,4,4],[4,4,5],
    [4,4,7],[4,4,9],[4,4,10],[4,5,5],[4,7,7],[4,7,9],
    [4,9,9],[4,10,10],[5,5,5],[7,7,7],[7,7,8],[7,7,9],
    [7,8,8],[7,9,9],[8,8,8],[8,8,9],[8,8,10],[8,9,9],
    [8,9,10],[8,10,10],[9,9,9],[9,9,10],[9,10,10],[10,10,10]]
    values = [4,1,1,1,1,1,
    1,1,1,1,1,6,
    1,1,1,2,1,1,
    1,1,1,2,7,2,
    1,1,1,2,1,1,
    1,1,4,4,1,1,
    1,1,6,1,1,1,
    1,1,4,1,1,3]
    lit_cp_tens["ex02_n10_d3_17Fan-Exa4.3"] = hcat(indices,values)

    #17Fan, Example 4.4
    indices = [[1,1,1,1],[1,1,1,10],[1,1,10,10],[1,10,10,10],[2,2,2,2],[2,2,2,4],
    [2,2,2,8],[2,2,2,9],[2,2,4,4],[2,2,4,8],[2,2,4,9],[2,2,8,8],
    [2,2,8,9],[2,2,9,9],[2,4,4,4],[2,4,4,8],[2,4,4,9],[2,4,8,8],
    [2,4,8,9],[2,4,9,9],[2,8,8,8],[2,8,8,9],[2,8,9,9],[2,9,9,9],
    [4,4,4,4],[4,4,4,8],[4,4,4,9],[4,4,8,8],[4,4,8,9],[4,4,9,9],
    [4,8,8,8],[4,8,8,9],[4,8,9,9],[4,9,9,9],[5,5,5,5],[5,5,5,7],
    [5,5,5,9],[5,5,7,7],[5,5,7,9],[5,5,9,9],[5,7,7,7],[5,7,7,9],
    [5,7,9,9],[5,9,9,9],[6,6,6,6],[6,6,6,7],[6,6,6,9],[6,6,6,10],
    [6,6,7,7],[6,6,7,9],[6,6,9,9],[6,6,10,10],[6,7,7,7],[6,7,7,9],
    [6,7,9,9],[6,9,9,9],[6,10,10,10],[7,7,7,7],[7,7,7,9],[7,7,9,9],
    [7,9,9,9],[8,8,8,8],[8,8,8,9],[8,8,8,10],[8,8,9,9],[8,8,9,10],
    [8,8,10,10],[8,9,9,9],[8,9,9,10],[8,9,10,10],[8,10,10,10],[9,9,9,9],
    [9,9,9,10],[9,9,10,10],[9,10,10,10],[10,10,10,10]]
    values = [1,1,1,1,6,2,
    2,2,2,1,1,2,
    1,2,2,1,1,1,
    1,1,2,1,1,2,
    6,2,2,2,1,2,
    2,1,1,2,2,1,
    1,1,1,1,1,1,
    1,1,3,1,1,1,
    1,1,1,1,1,1,
    1,1,1,4,2,2,
    2,8,3,1,3,1,
    1,3,1,1,1,12,
    1,1,1,4]
    lit_cp_tens["ex03_n10_d4_17Fan-Exa4.4"] = hcat(indices,values)

    #14Qi, Example 2(1) # In the article of 14Qi, they miss [2,10,10]->1
    indices = [[1,1,1],[1,1,5],[1,5,5],
    [2,2,2],[2,2,3],[2,2,6],[2,2,8],[2,2,9],[2,2,10],[2,3,3],[2,6,6],[2,6,9],[2,8,8],[2,8,10],[2,9,9],[2,10,10],
    [3,3,3],[3,3,4],[3,3,5],[3,4,4],[3,4,5],[3,5,5],
    [4,4,4],[4,4,5],[4,5,5],
    [5,5,5],[5,5,9],[5,9,9],
    [6,6,6],[6,6,9],[6,9,9],
    [7,7,7],[7,7,9],[7,7,10],[7,9,9],[7,9,10],[7,10,10],
    [8,8,8],[8,8,10],[8,10,10],
    [9,9,9],[9,9,10],[9,10,10],
    [10,10,10]]
    values = [1,1,1,
    5,1,1,1,1,1,1,1,1,1,1,1,1,
    3,1,1,1,1,1,
    2,1,1,
    4,1,1,
    2,1,1,
    2,1,1,1,1,1,
    2,1,1,
    5,1,1,
    4]
    lit_cp_tens["ex04_n10_d3_14Qi-Exa2-1"] = hcat(indices,values)

    #14Qi, Example 2(2)
    indices = [[1,1,1],[1,1,5],[1,1,10],[1,5,5],[1,5,10],[1,10,10],
    [2,2,2],[2,2,3],[2,2,8],[2,2,9],[2,2,10],[2,3,3],[2,3,9],[2,8,8],[2,9,9],[2,9,10],[2,10,10],
    [3,3,3],[3,3,4],[3,3,8],[3,3,9],[3,4,4],[3,4,8],[3,8,8],[3,8,9],[3,9,9],
    [4,4,4],[4,4,8],[4,8,8],
    [5,5,5],[5,5,8],[5,5,10],[5,8,8],[5,10,10],
    [8,8,8],[8,8,9],[8,9,9],
    [9,9,9],[9,9,10],[9,10,10],
    [10,10,10]]
    values = [2,1,1,1,1,1,
    5,1,1,2,1,1,1,1,2,1,1,
    6,1,2,2,1,1,2,1,2,
    2,1,1,
    3,1,1,1,1,
    6,1,1,
    6,1,1,
    4]
    lit_cp_tens["ex05_n10_d3_14Qi-Exa2-2"] = hcat(indices,values)

    #14Qi, Example 3(2)
    indices = [[1,1,1,1],[1,1,1,2],[1,1,1,3],[1,1,1,5],[1,1,1,7],[1,1,1,8],[1,1,1,9],[1,1,2,2],[1,1,2,7],
    [1,1,3,3],[1,1,3,8],[1,1,3,9],[1,1,5,5],[1,1,7,7],[1,1,8,8],[1,1,8,9],[1,1,9,9],[1,2,2,2],[1,2,2,7],[1,2,7,7],[1,3,3,3],
    [1,3,3,8],[1,3,3,9],[1,3,8,8],[1,3,8,9],[1,3,9,9],[1,5,5,5],[1,7,7,7],[1,8,8,8],[1,8,8,9],[1,8,9,9],[1,9,9,9],
    [2,2,2,2],[2,2,2,3],[2,2,2,6],[2,2,2,7],[2,2,3,3],[2,2,3,6],[2,2,6,6],[2,2,6,7],[2,2,7,7],[2,3,3,3],[2,3,3,6],[2,3,6,6],[2,6,6,6],[2,6,6,7],[2,6,7,7],[2,7,7,7],
    [3,3,3,3],[3,3,3,6],[3,3,3,8],[3,3,3,9],[3,3,6,6],[3,3,8,8],[3,3,8,9],[3,3,9,9],[3,6,6,6],[3,8,8,8],[3,8,8,9],[3,8,9,9],[3,9,9,9],
    [4,4,4,4],[4,4,4,9],[4,4,9,9],[4,9,9,9],
    [5,5,5,5],
    [6,6,6,6],[6,6,6,7],[6,6,7,7],[6,7,7,7],
    [7,7,7,7],[7,7,7,9],[7,7,7,10],[7,7,9,9],[7,7,9,10],[7,7,10,10],[7,9,9,9],[7,9,9,10],[7,9,10,10],[7,10,10,10],
    [8,8,8,8],[8,8,8,9],[8,8,9,9],[8,9,9,9],
    [9,9,9,9],[9,9,9,10],[9,9,10,10],[9,10,10,10],
    [10,10,10,10]]
    values = [9,1,2,1,1,2,2,1,1,
    2,1,1,1,1,2,1,2,1,1,1,2,
    1,1,1,1,1,1,1,2,1,1,2,
    6,1,2,2,1,1,2,1,2,1,1,1,2,1,1,2,
    8,1,2,2,1,2,1,2,1,2,1,1,2,
    1,1,1,1,
    1,
    4,1,1,1,
    6,1,1,1,1,1,1,1,1,1,
    6,2,2,2,
    9,1,1,1,
    2]
    lit_cp_tens["ex06_n10_d4_14Qi-Exa3-2"] = hcat(indices,values)
    
    #14Qi, Example 3(3)
    indices = [[1,1,1,1],[1,1,1,5],[1,1,1,6],[1,1,1,8],[1,1,1,9],[1,1,1,10],
    [1,1,5,5],[1,1,5,6],[1,1,5,8],[1,1,5,9],[1,1,5,10],[1,1,6,6],[1,1,6,8],[1,1,6,9],[1,1,8,8],[1,1,9,9],[1,1,9,10],[1,1,10,10],
    [1,5,5,5],[1,5,5,6],[1,5,5,8],[1,5,5,9],[1,5,5,10],[1,5,6,6],[1,5,6,8],[1,5,6,9],[1,5,8,8],[1,5,9,9],[1,5,9,10],[1,5,10,10],
    [1,6,6,6],[1,6,6,8],[1,6,6,9],[1,6,8,8],[1,6,9,9],[1,8,8,8],[1,9,9,9],[1,9,9,10],[1,9,10,10],[1,10,10,10],
    [2,2,2,2],[2,2,2,5],[2,2,2,6],[2,2,2,9],[2,2,5,5],[2,2,5,6],[2,2,5,9],[2,2,6,6],[2,2,6,9],[2,2,9,9],
    [2,5,5,5],[2,5,5,6],[2,5,5,9],[2,5,6,6],[2,5,6,9],[2,5,9,9],[2,6,6,6],[2,6,6,9],[2,6,9,9],[2,9,9,9],
    [3,3,3,3],[3,3,3,4],[3,3,3,9],[3,3,3,10],[3,3,4,4],[3,3,4,9],[3,3,9,9],[3,3,9,10],[3,3,10,10],
    [3,4,4,4],[3,4,4,9],[3,4,9,9],[3,9,9,9],[3,9,9,10],[3,9,10,10],[3,10,10,10],
    [4,4,4,4],[4,4,4,9],[4,4,9,9],[4,9,9,9],
    [5,5,5,5],[5,5,5,6],[5,5,5,8],[5,5,5,9],[5,5,5,10],[5,5,6,6],[5,5,6,8],[5,5,6,9],[5,5,8,8],[5,5,8,9],[5,5,9,9],[5,5,9,10],[5,5,10,10],
    [5,6,6,6],[5,6,6,8],[5,6,6,9],[5,6,8,8],[5,6,9,9],[5,8,8,8],[5,8,8,9],[5,8,9,9],[5,9,9,9],[5,9,9,10],[5,9,10,10],[5,10,10,10],
    [6,6,6,6],[6,6,6,8],[6,6,6,9],[6,6,8,8],[6,6,9,9],[6,8,8,8],[6,9,9,9],
    [8,8,8,8],[8,8,8,9],[8,8,9,9],[8,9,9,9],
    [9,9,9,9],[9,9,9,10],[9,9,10,10],[9,10,10,10],
    [10,10,10,10]]
    values = [18,6,4,2,4,2,
    6,2,1,2,1,4,1,1,2,4,1,2,
    6,2,1,2,1,2,1,1,1,2,1,1,
    4,1,1,1,1,2,4,1,1,2,
    6,2,2,2,2,1,1,2,1,2,
    2,1,1,1,1,1,2,1,1,2,
    4,1,2,1,1,1,2,1,1,
    1,1,1,2,1,1,1,
    2,1,1,1,
    26,6,3,7,2,6,1,2,3,1,7,1,2,
    6,1,2,1,2,3,1,1,7,1,1,2,
    18,2,4,2,4,2,4,
    8,1,1,1,
    24,3,3,3,
    8]
    lit_cp_tens["ex07_n10_d4_14Qi-Exa3-3"] = hcat(indices,values)
    
    return lit_cp_tens
end
function generate_lit_cp_tens(save_dir::String)
    !isdir(save_dir) ? mkdir(save_dir) : 0
    lit_cp_tens = generate_lit_cp_tens()
    for k in keys(lit_cp_tens)
        mat = lit_cp_tens[k]
        matrix_IO.save_tms(mat,save_dir*k*"_.csv")
    end
end

### Literature non Examples
function generate_lit_non_cp_tens()
    lit_non_cp_tens = Dict()
    #17Fan, Example 4.1
    indices = [[1,1,1],[1,1,2],[1,2,5],[1,2,9],[1,3,8],[1,5,11],
    [2,2,2],[2,2,3],[2,3,5],[2,3,7],[2,4,6],[2,4,11],
    [3,3,3],[3,3,6],[3,3,9],[3,4,6],[3,5,7],[3,5,10],
    [3,7,10],[3,8,9],[3,9,11],[3,10,10],[4,4,4],[4,4,5],
    [4,4,9],[4,4,11],[4,5,7],[4,5,9],[4,6,10],[4,7,9],
    [4,9,9],[4,10,11],[5,5,5],[5,5,10],[5,7,9],[5,8,8],
    [5,9,10],[5,10,10],[5,11,11],[6,6,6],[7,7,7],[7,7,9],
    [7,8,10],[7,9,11],[8,8,8],[8,8,10],[8,8,11],[8,9,10],
    [8,9,11],[8,10,10],[9,9,9],[9,9,10],[9,10,10],[10,10,10],
    [10,11,11],[11,11,11]]
    values = [5,3,1,3,2,1,
    7,2,1,1,1,1,
    6,1,2,1,2,1,
    1,3,1,1,5,2,
    2,1,2,2,1,1,
    1,1,5,1,1,1,
    1,1,1,4,3,1,
    2,2,6,1,1,2,
    3,1,4,1,1,3,1,5]
    lit_non_cp_tens["ex01_n11_d3_17Fan-Exa4.1"] = hcat(indices,values)

    #17Fan, Example 4.2
    indices = [[1,1,1,1,1],[1,1,2,2,2],[1,2,4,4,6],[1,5,6,6,8],[1,6,7,8,8],
    [2,2,2,2,2],[2,2,3,5,7],[2,3,4,6,8],[2,5,6,7,7],[2,6,7,7,8],
    [3,3,3,3,3],[3,4,5,6,6],[3,5,7,8,8],[3,6,6,7,8],[3,7,7,7,7],
    [4,4,4,4,4],[4,5,5,6,8],[4,6,6,7,8],[4,6,7,8,8],[4,7,7,7,7],
    [5,5,5,5,5],[5,5,6,6,7],[5,6,7,7,8],[5,7,7,8,8],[5,7,8,8,8],
    [6,6,6,6,6],[6,6,6,6,7],[6,7,8,8,8],[6,8,8,8,8],[7,7,7,7,7],
    [7,7,7,8,8],[7,7,8,8,8],[8,8,8,8,8]]
    values = [3,1,1,1,2,
    5,2,1,1,1,
    3,1,2,1,2,
    6,2,2,1,2,
    5,1,1,1,1,
    6,2,2,1,7,
    2,1,8]
    lit_non_cp_tens["ex02_n8_d5_17Fan-Exa4.2"] = hcat(indices,values)
    return lit_non_cp_tens
end
function generate_lit_non_cp_tens(save_dir::String)
    !isdir(save_dir) ? mkdir(save_dir) : 0
    lit_non_cp_tens = generate_lit_non_cp_tens()
    for k in keys(lit_non_cp_tens)
        mat = lit_non_cp_tens[k]
        matrix_IO.save_tms(mat,save_dir*k*"_.csv")
    end
end

### Random Example
function generate_random_cp_tens()
    n = [6,8,10,12,14] # various dimensions of tensors
    m = [4,6,8] # various orders of tensors
    p = [0.4,0.6,0.8,0.9,0.98] # nzd: the sparsity of tensors
    repeattime = 5; # the repeat time for per dimension, order and sparsity

    N = repeat(n,inner = length(m)*length(p)) # N = [6 repeats 15 times, 8 repeats 15 times,..., 14 repeats 15 times]
    N = repeat(N,outer=repeattime) #[6 repeats 15 times, 8 repeats 15 times,...,14 repeats 15 times, 6 repeats 15 times, 8 repeats 15 times,...,14 repeats 15 times, ...]
    M = repeat(m,inner=length(p)) #[4 repeats 5 times, 6 repeats 5 times, 8 repeat 5 times]
    M = repeat(M,outer=length(n)*repeattime) #[4 repeats 5 times, 6 repeats 5 times, 8 repeat 5 times, ...]
    P = repeat(p,outer=length(n)*length(m)*repeattime) #[the whole p repeats 75 times]
    #N,M,P:(6,4,0.4),...,(6,4,0.98),(6,6,0.4),...,(6,6,0.98),(6,8,0.4),...,(6,8,0.98),(8,4,0.4) and then all repeat repeattimes
    return generate_random_cp_tens(N,M,P)
end
function generate_random_cp_tens(save_dir::String)
    !isdir(save_dir) ? mkdir(save_dir) : 0
    rand_cp_tms = generate_random_cp_tens()
    for k in keys(rand_cp_tms)
        tms = rand_cp_tms[k]
        matrix_IO.save_tms(tms,save_dir*k*"_.csv")
    end
end
function generate_random_cp_tens(N::Vector{Int64}, M::Vector{Int64}, P)
    rand_cp_tens = Dict()
    count = 0
    for (n,m,p) ∈ [zip(N,M,P)...]
        count += 1
        tms = generate_random_cp_tens((n,m),p)
        indices = moments.get_subscripts(m,n)
        values = [tms[item] for item in indices]
        tms = hcat(indices,values)
        
        ext = lpad(count,3,'0') # the number of test examples, pad '0' at the left of the count to obtain a 3 length string
        tms_name = "R_ex$(ext)_n$(n)_d$(m)_p$(p)_"
        rand_cp_tens[tms_name] = tms
    end
    return rand_cp_tens
end
function generate_random_cp_tens(N::Vector{Int64}, M::Vector{Int64}, P, save_dir::String)
    !isdir(save_dir) ? mkdir(save_dir) : 0
    rand_cp_tms = generate_random_cp_tens(N,M,P)
    for k in keys(rand_cp_tms)
        tms = rand_cp_tms[k]
        matrix_IO.save_tms(tms,save_dir*k*"_.csv")
    end
end
function generate_random_cp_tens((n,m)::Tuple{Int64,Int64},p)
    # generate the tms corresponding to a symmetric binary tensor, the number of nonzeros is (p*(binomial(n+m-1,m)-n))+n
    M = gen_random_supp_graph(n,m,p)

    """
    # Do not use the code below, which has not work properly yet, we wish it generate the cp decomposition according to the support multi-hypergraph
    mc = moments.get_maximal_cliques(M,n,m)
    nums = length(mc)
    a_s = [embed(round.(rand(length(mc[k])),digits=2),mc[k],n) for k in 1:nums]
    a_s = [a_s..., [m,p,repeat([0],n-2)...]] # the lost column record the order m and the nzd p of tensors
    return a_s
    """
    return M
end
embed(v,β,n) = [i ∈ β ? popfirst!(v) : 0.0 for i in 1:n]

function gen_random_supp_graph(n, m, p)
    # For almost all the examples, it hard to generate a cp tensor randomly
    subscripts = moments.get_subscripts(m, n)
    # Determin which position should be zero randomly
    n1 = ceil(Int64,p*(length(subscripts)-n)) # the diagonal entry of our random generated tensor will always be 1

    M = Dict()
    rand_0_1_vec = shuffle([ones(Int64, n1)...,zeros(Int64,length(subscripts)-n1-n)...]) #shuffle return a random arrangement of the elements in a list
    for i in 1:length(subscripts)
        if length(Set(subscripts[i])) !=1
            M[subscripts[i]] = popfirst!(rand_0_1_vec)
        else
            M[subscripts[i]] = 1
        end
    end
    #mc = moments.get_maximal_cliques(M,n,m)
    #check = all(a->!isempty(a),[filter(c->issubset(item,c),mc) for item in moments.get_nonzero_entries(M)])
    
    return M

end

function gen_random_cp_supp_graph(n, m, p)
    # Dot not use this function, not work yet, it hard to generate a cp tensor randomly，while loop cannot terminate
        subscripts = moments.get_subscripts(m, n)
    # Determin which position should be zero randomly
    n1 = ceil(Int64,p*(length(subscripts)-n)) # the diagonal entry of our random generated tensor will always be 1
    check = false
    M = Dict()
    while !check
        M = Dict()
        rand_0_1_vec = shuffle([ones(Int64, n1)...,zeros(Int64,length(subscripts)-n1-n)...])
        for i in 1:length(subscripts)
            if length(Set(subscripts[i])) !=1
                M[subscripts[i]] = popfirst!(rand_0_1_vec)
            else
                M[subscripts[i]] = 1
            end
        end
        mc = moments.get_maximal_cliques(M,n,m)
        check = all(a->!isempty(a),[filter(c->issubset(item,c),mc) for item in moments.get_nonzero_entries(M)])
    end
    return M

end

end