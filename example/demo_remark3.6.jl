#----------------------- Initialization ---------------------
src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"cp_model.jl")
include(src_dir*"extract_atoms.jl")
include(src_dir*"my_extract_atoms.jl")

#Completely positive tensor with [[1,2,3],[1,4]] cliques
#But when genereted, we use v1=[xx,xx,0,0],v2=[xx,0,xx,0],v3=[0,xx,xx,0],v4=[xx,0,0,xx]
#and the resulted decomposition by the algorithm may be v1=[0,xx,xx,0],[xx,xx,xx,0],[xx,0,xx,0];[xx,0,0,xx]
using Random
Random.seed!(1234); # generate repeatable random numbers

d = 2; n = 4;
v1 = rand(4)
v1[3] = 0.0; v1[4] = 0.0
v2 = rand(4)
v2[2] = 0.0; v2[4] = 0.0
v3 = rand(4)
v3[1] = 0.0; v3[4] = 0.0
v4 = rand(4)
v4[2] = 0.0; v4[3] = 0.0
v5 = zeros(4); v5[1] = 2 # the last column records the order of this tensor
vmat = hcat(v1,v2,v3,v4,v5)

A = matrix_IO.get_tms_from_factors(vmat) # the tms corresponding to the tensor

t₀ = Int(ceil( (d+1)/2 )) # the initial relaxation order
len = binomial(n+t₀,t₀)
G = rand(len,len)
F = G'*G

lvt = t₀+2

#-- ideal dense model----------------------------------------------------------------------
#=ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"id")

n_id, t_id, mom_vals_id = extract_atoms.proc_mom(mom_mat_id)

try
    println("----id extract----")
    ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)

    cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
    A_ext_id                = extract_atoms.recon_ten_id(cents_id, weights_id,d)

    error_id = sum([abs(A_ext_id[key]-A[key]) for key in keys(A)])
    println("id_recover_error: ", error_id)
catch
    println("No atom could be extracted")
end=#

#--ideal sparse model-----------------------------------------------------------------------
cpu = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end

n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(mom_mat_sp)

#=try
    println("----sp extract-----")
    ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-4, true)

    cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
    println("sp_recover_error: ", error_sp)
catch
    println("No atom could be extracted")
end=#

#try
    println("----sp extract-----")
    ext_atoms_sp            = my_extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-6)

    cents_sp, weights_sp    = my_extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = my_extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    err = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
    println("sp_recover_error: ", err)
#catch
    #println("No atom could be extracted")
#end