module tensorGraphs

function maximal_cliques(A,n,d)
    # Obtain all the subscripts of zero entries
    zero_subscripts = sort([item for item in keys(A) if A[item] ==0])

    cand = [[1:n...]]
    for item in zero_subscripts
        # spilt procedure, the nodes in item cannot beglone to the same clique
        add_cand = []
        for s in cand
            if issubset(item,s)
                for k in Set(item)
                    temp = s[:]
                    s_new = filter!(a->a!=k,temp)
                    push!(add_cand,s_new)
                end
            end
        end
        # Delete the candidate set that cannot be a clique from cand
        filter!(a->!issubset(item,a), cand)
        # Add the new possible candidate set that may be a clique to cand, and ensure that the sets in cand do not contain eachother
        while !isempty(add_cand)
            temp = pop!(add_cand)
            check = true
            for s in cand
                if issubset(temp,s)
                    check=false
                    break
                end
            end
            if check
                push!(cand,temp)
            end
        end
    end
    return cand
end

end