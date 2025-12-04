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
#len = binomial(t0+n,t0)
#G = rand(len,len)
#F = G'*G 
F_list = readdir(assets_dir*"results\\literature_non\\t$(t0)\\F\\")
F = Matrix(matrix_IO.load_mat(assets_dir*"results\\literature_non\\t$(t0)\\F\\"*F_list[1]))

lvt = t0

#-- ideal dense model----------------------------------------------------------------------
total_time = @elapsed begin
cpu = @elapsed begin ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"id") end
#@time begin ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"id") end

n_id, t_id, mom_vals_id = extract_atoms.proc_mom(mom_mat_id)
check = my_extract_atoms.rank_check(n_id, t_id, d, mom_vals_id, 1e-4)
try
    println("----id extract----")
    ext_attimeoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)

    cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
    A_ext_id                = extract_atoms.recon_ten_id(cents_id, weights_id,d)

    error_id = sum([abs(A_ext_id[key]-A[key]) for key in keys(A)])
    println("id_recover_error: ", error_id)
catch
    println("No atom could be extracted")
end
end