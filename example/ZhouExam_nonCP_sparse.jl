src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"matrix_IO.jl")
include(src_dir*"cp_model.jl")
include(src_dir*"extract_atoms.jl")
include(src_dir*"my_extract_atoms.jl")

#-----------------------------------Load Data------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
mat_types = ["literature", "literature_non", "random"]

mat_dir = data_dir*mat_types[2]*"\\"
data_matrix_dict = matrix_IO.load_tmses(mat_dir)
dmdk = data_matrix_dict.keys
dmdv = data_matrix_dict.vals

k=2
A_name = dmdk[k]
A = dmdv[k]

#----------------------------------SDP relaxation--------------------------------------
n = parse(Int64,split(A_name,"_")[2][2:end])
d = parse(Int64,split(A_name,"_")[3][2:end])
t0 = Int(ceil((d+1)/2))
len = binomial(t0+n,t0)
G = rand(len,len)
F = G'*G 

lvt = t0 

#-- ideal sparse model----------------------------------------------------------------------
#total_time = @elapsed begin
@elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end
cpu = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end

n_sp, t_sp, mom_vals_sp = my_extract_atoms.proc_mom(mom_mat_sp)
tol = 1e-6
check = my_extract_atoms.rank_check(n_sp, t_sp, d, mom_vals_sp, tol)
try
    println("----sp extract-----")
    ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, tol)

    cents_sp, weights_sp    = my_extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = my_extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
    println("sp_recover_error: ", error_sp)
catch
    println("No atom could be extracted")
end
#end