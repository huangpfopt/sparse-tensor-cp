src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"matrix_IO.jl")
include(src_dir*"cp_model.jl")
include(src_dir*"extract_atoms.jl")
include(src_dir*"my_extract_atoms.jl")

using Random
Random.seed!(1234);
#-----------------------------------Load Data------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
mat_types = ["literature", "literature_non", "random"]

mat_dir = data_dir*mat_types[1]*"\\"
data_matrix_dict = matrix_IO.load_tmses(mat_dir)
dmdk = data_matrix_dict.keys
dmdv = data_matrix_dict.vals

k=7 # for dense model, we can only compute the case of d=3, when d=4, it occurs Mosek.MosekError(1051, ""), that is out of space
A_name = dmdk[k]
A = dmdv[k]


#----------------------------------SDP relaxation--------------------------------------
n = parse(Int64,split(A_name,"_")[2][2:end])
d = parse(Int64,split(A_name,"_")[3][2:end])
t0 = Int(ceil((d+1)/2))
#len = binomial(t0+n,t0)
#G = randn(len,len)
#F = G'*G 
#F = F/sqrt(sum(F.*F))
F_list = readdir(assets_dir*"results\\literature\\t$(t0)\\F\\")
F = Matrix(matrix_IO.load_mat(assets_dir*"results\\literature\\t$(t0)\\F\\"*F_list[3]))

lvt = t0

#-- ideal dense model----------------------------------------------------------------------
cpu = @elapsed begin ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"id") end

tol = 1e-6
n_id, t_id, mom_vals_id = extract_atoms.proc_mom(mom_mat_id)
check = my_extract_atoms.rank_check(n_id, t_id, d, mom_vals_id, tol)
try
    println("----id extract----")
    ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, tol, true)

    cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
    A_ext_id                = extract_atoms.recon_ten_id(cents_id, weights_id,d)

    error_id = sum([abs(A_ext_id[key]-A[key]) for key in keys(A)])
    println("id_recover_error: ", error_id)
catch
    println("No atom could be extracted")
end