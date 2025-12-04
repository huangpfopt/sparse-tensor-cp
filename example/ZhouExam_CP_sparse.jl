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
data_tensor_dict_cp = matrix_IO.load_tmses(mat_dir)
dmdk = data_tensor_dict_cp.keys
dmdv = data_tensor_dict_cp.vals

k=7
A_name = dmdk[k]
A = dmdv[k]
#A_value = dmdv[k]
#mat = hcat(eval.(Meta.parse.(A_value.tms_inds)), A_value.tms_vals)
#A = matrix_IO.get_tms_from_indices(mat)

#----------------------------------SDP relaxation--------------------------------------
n = parse(Int64,split(A_name,"_")[2][2:end])
d = parse(Int64,split(A_name,"_")[3][2:end])
t0 = Int(ceil((d+1)/2))
#len = binomial(t0+n,t0)
#G = rand(len,len)
#F = G'*G 
F_list = readdir(assets_dir*"results\\literature\\t$(t0)\\F\\")
F = Matrix(matrix_IO.load_mat(assets_dir*"results\\literature\\t$(t0)\\F\\"*F_list[3]))

#lvt = t0+1 # for k=3,6,7 lvt+1 can obtain higher recover accuracy(by hom or robust extract methods)
lvt = t0

#-- ideal sparse model----------------------------------------------------------------------
#total_time = @elapsed begin
@elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end
cpu_sp = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end
println("cpu_sp: ", cpu_sp)

n_sp, t_sp, mom_vals_sp = my_extract_atoms.proc_mom(mom_mat_sp)
tol = 1e-6
check = my_extract_atoms.rank_check(n_sp, t_sp, d, mom_vals_sp, tol)

#=try
    println("----sp extract-----")
    ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, tol)

    cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)]) # Since random numbers are used, there may be slight variations
    println("sp_recover_error: ", error_sp)
catch
    println("No atom could be extracted")
end=#

try
    println("----sp extract-----")
    #ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-6)
    #ts = [Int(ceil(d/2))-1+length(rs[2])-1 for rs in check[2]]
    #ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, ts, mom_vals_sp, 1e-6)
    ext_atoms_sp            = my_extract_atoms.get_atoms_robust(n_sp, t_sp, mom_vals_sp, 1e-6)

    cents_sp, weights_sp    = my_extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = my_extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    err = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
    println("sp_recover_error: ", err)
catch
    println("No atom could be extracted")
end

#end