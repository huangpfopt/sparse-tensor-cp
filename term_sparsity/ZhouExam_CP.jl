using TSSOS
using DynamicPolynomials
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

mat_dir = data_dir*mat_types[1]*"\\"
data_tensor_dict_cp = matrix_IO.load_tmses(mat_dir)
dmdk = data_tensor_dict_cp.keys
dmdv = data_tensor_dict_cp.vals


k=7 # for all cases, the sdp solved time of sp model always less than of tssos, however if included the modeling time, only from d=4, ideal sparse show obvious advantage
# as solving the primal problem, for k=2,4,5, lvt = 2, tssos in general cannot obtain correct extract atoms, it need to increase lvt to 3
A_name = dmdk[k]
A = dmdv[k]

#----------------------------------SDP relaxation--------------------------------------
n = parse(Int64,split(A_name,"_")[2][2:end])
d = parse(Int64,split(A_name,"_")[3][2:end])
t0 = Int(ceil((d+1)/2))
#len = binomial(t0+n,t0)
#G = rand(len,len)
#F = G'*G 
F_list = readdir(assets_dir*"results\\literature\\t$(t0)\\F\\")
F = Matrix(matrix_IO.load_mat(assets_dir*"results\\literature\\t$(t0)\\F\\"*F_list[3]))

lvt = t0 # begin from lvt=4, TSSOS may also have LoadError: Mosek.MosekError(1051, "") or take very long time,  while ideal sparse model still works 

tol = 1e-6
#=#----------------------------ideal sparse model---------------------------------------------
total_sp = @elapsed begin
cpu_sp = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end
println("cpu_sp: ", cpu_sp)

n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(mom_mat_sp)

check = my_extract_atoms.rank_check(n_sp, t_sp, d, mom_vals_sp, tol)

println("----sp extract-----")
# extract_atoms may sometimes obtain higher recover accuracy than my_extract_atoms
ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, tol)

cents_sp, weights_sp    = my_extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)

A_ext_sp                = my_extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
println("sp_recover_error: ", error_sp)
end

# show(ξₜᶜᵖˢᵖ) # show(model) will print a summary of the problem=#
#----------------TSSOS model--------------------------------------
#total_tssos = @elapsed begin

@elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F) end
cpu_tssos = @elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F) end
println("cpu_tssos: ", cpu_tssos)
mom_vals_tssos = mom_mat_tssos
MomMat = TSSOS.get_moment_matrix(mom_vals_tssos, info1)

println("----tssos extract-----")
sol = extract_solutions_robust(n, lvt, MomMat[1]; tol=tol)
#@polyvar x[1:n]
#sol = extract_solutions([zeros(n)'*x],x,lvt, 0, MomMat[1]; tol=tol) # when k=2, no atoms may cannot be extracted
# compute the weights
I₂ₜ = moments.make_mon_expo(n, d; isle=false)
mom = zeros(length(I₂ₜ), length(sol))
for i in 1:length(sol) 
    mom[:,i] = map(a->prod(sol[i].^a),I₂ₜ)
end
weights = mom \ [A[item] for item in moments.get_subscripts(d,n)]
A_ext_tssos = my_extract_atoms.recon_ten_id(sol,weights,d)
error_tssos = sum([abs(A_ext_tssos[key]-A[key]) for key in keys(A)])
println("tssos_recover_error: ", error_tssos)
#end
