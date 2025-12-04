using TSSOS
using Random
Random.seed!(1234);
#----------------------- Initialization ---------------------
src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"extract_atoms.jl")
include(src_dir*"my_extract_atoms.jl")
include(src_dir*"cp_model.jl")

tssos_dir = dirname(@__FILE__)*"\\"
include(tssos_dir*"tssos_model.jl")

#-----------------------------------Load Data------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
mat_types = ["literature", "literature_non", "random"]

mat_dir = data_dir*mat_types[2]*"\\"
data_tensor_dict_noncp = matrix_IO.load_tmses(mat_dir)
dmdk = data_tensor_dict_noncp.keys
dmdv = data_tensor_dict_noncp.vals

k=2 # simial to the cp case, only when k=2 (that is lvt≥3), sp has obvious advantage. As for sdp solve time, sp always better
A_name = dmdk[k]
A = dmdv[k]

#----------------------------------SDP relaxation--------------------------------------
n = parse(Int64,split(A_name,"_")[2][2:end])
d = parse(Int64,split(A_name,"_")[3][2:end])
t0 = Int(ceil((d+1)/2))
#len = binomial(t0+n,t0)
#G = rand(len,len)
#F = G'*G 
F_list = readdir(assets_dir*"results\\literature_non\\t$(t0)\\F\\")
F = Matrix(matrix_IO.load_mat(assets_dir*"results\\literature_non\\t$(t0)\\F\\"*F_list[1]))


lvt = t0

#----------------TSSOS model--------------------------------------
@elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F) end
cpu_tssos = @elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F) end
println("cpu_tssos: ", cpu_tssos)
mom_vals_tssos = mom_mat_tssos
MomMat = TSSOS.get_moment_matrix(mom_vals_tssos, info1)

tol = 1e-6
sol = extract_solutions_robust(n, lvt, MomMat[1]; tol=tol)


#=#----------------------------ideal sparse model---------------------------------------------
cpu_sp = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end

n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(mom_mat_sp)

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
end=#