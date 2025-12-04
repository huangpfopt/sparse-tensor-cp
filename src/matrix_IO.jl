module matrix_IO
using CSV, DataFrames

proj_dir = dirname(@__FILE__)*"\\" #return the absolute path of current file
include(proj_dir*"moments.jl")

### Matrices
save_mat(M, save_path::String) = CSV.write(save_path, DataFrame(M, :auto))
load_mat(load_path::String) = CSV.read(load_path, DataFrame)
function load_mats(load_dir::String)
    if load_dir[end-3:end] == ".csv"
        return load_mat(load_dir)
    else
        list_o_mats = [s for s in readdir(load_dir) if (contains(s,".csv") && contains(s,"ex"))]
        list_o_mat_names = [m[1:end-4] for m in list_o_mats]
        mats = [load_mat(load_dir*mat) for mat in list_o_mats]
        return sort(Dict(zip(list_o_mat_names, mats))) #return an orderedDict, ordered by keys
    end
end

### tms
save_tms(M, save_path::String) = CSV.write(save_path, DataFrame(tms_inds=M[:,1], tms_vals=M[:,2]))
function load_tms(load_path::String)
    mat = CSV.read(load_path, DataFrame)
    mat = hcat(eval.(Meta.parse.(mat.tms_inds)), mat.tms_vals)
    A = matrix_IO.get_tms_from_indices(mat)
    return A
end
"""Loads all the tensors data stored as .csv's from a directory"""
function load_tmses(load_dir::String)
    if load_dir[end-3:end] == ".csv"
        return load_tms(load_dir)
    else
        list_o_mats = [s for s in readdir(load_dir) if (contains(s,".csv") && contains(s,"ex"))]
        list_o_ten_names = [m[1:end-4] for m in list_o_mats]
        tens = [load_tms(load_dir*mat) for mat in list_o_mats]
        return sort(Dict(zip(list_o_ten_names, tens)))
    end
end

"""从factors中获取张量对应的tms"""
function get_tms_from_factors(mats)
    n = size(mats,1)
    d = Int64(mats[1,end]) # the first element of the last column in mats records the order of tensors

    A = Dict()
    subscripts = moments.get_subscripts(d,n)
    for item in subscripts
        A[item] = sum([prod([mats[i,j] for i in item]) for j in 1:size(mats,2)-1])
    end
    return A
end

"""从csv中存储的indices和values里恢复出tms"""
function get_tms_from_indices(mats)
    indices = mats[:,1]
    values = mats[:,2].+0.0
    n = max([max(item...) for item in indices]...)
    d = length(indices[1])
    A = Dict(zip(indices,values))
    subscripts = moments.get_subscripts(d,n)
    for item in subscripts
        if !(item in keys(A))
            A[item] = 0.0
        end
    end
    return A
end

### moments
function save_moments(ext_mom, save_path::String)
    df = DataFrame(mom_inds = ext_mom[:,1], mom_vals = ext_mom[:,2])
    CSV.write(save_path, df)
end
function load_moments(load_path::String)
    mat = CSV.read(load_path, DataFrame)
    mat = hcat(eval.(Meta.parse.(mat.mom_inds)), mat.mom_vals)
    return mat   
end

end