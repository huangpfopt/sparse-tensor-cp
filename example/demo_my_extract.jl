#----------------------- Initialization ---------------------
src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"cp_model.jl")
include(src_dir*"my_extract_atoms.jl")

#Completely positive tensor with [[1,2,3],[2,3,4],[1,4,5]] cliques
using Random
Random.seed!(1234); # generate repeatable random numbers

d = 3; n = 5;
v1 = rand(5)
v1[4] = 0.0; v1[5] = 0.0
v2 = rand(5)
v2[1] = 0.0; v2[5] = 0.0
v3 = rand(5)
v3[2] = 0.0; v3[3] = 0.0
v4 = rand(5)
v4[2] = 0.0; v4[3] = 0.0
v5 = zeros(5); v5[1] = 3 # the last column records the order of this tensor
vmat = hcat(v1,v2,v3,v4,v5)

A = matrix_IO.get_tms_from_factors(vmat) # the tms corresponding to the tensor

t₀ = Int(ceil( (d+1)/2 )) # the initial relaxation order
len = binomial(n+t₀,t₀)
G = rand(len,len)
F = G'*G

lvt = t₀

#=#-- ideal dense model----------------------------------------------------------------------
ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"id")

n_id, t_id, mom_vals_id = my_extract_atoms.proc_mom(mom_mat_id)

check = my_extract_atoms.rank_check(n_id, t_id, d, mom_vals_id, 1e-4)
try
    println("----id extract----")
    ext_atoms_id            = my_extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4)

    cents_id, weights_id    = my_extract_atoms.ext_centers_weights(ext_atoms_id)
    A_ext_id                = my_extract_atoms.recon_ten_id(cents_id, weights_id,d)

    err = sum([abs(A_ext_id[key]-A[key]) for key in keys(A)])
    println("id_recover_error: ", err)
catch
    println("No atom could be extracted")
end=#

#--ideal sparse model-----------------------------------------------------------------------
cpu = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end

n_sp, t_sp, mom_vals_sp = my_extract_atoms.proc_mom(mom_mat_sp)

check = my_extract_atoms.rank_check(n_sp, t_sp, d, mom_vals_sp, 1e-4)
try
    println("----sp extract-----")
    ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-4)

    cents_sp, weights_sp    = my_extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = my_extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    err = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
    println("sp_recover_error: ", err)
catch
    println("No atom could be extracted")
end