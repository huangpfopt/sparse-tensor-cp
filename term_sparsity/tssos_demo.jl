using TSSOS
using DynamicPolynomials
#----------------------- Initialization ---------------------
src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"extract_atoms.jl")
include(src_dir*"cp_model.jl")
include(src_dir*"my_extract_atoms.jl")

tssos_dir = dirname(@__FILE__)*"\\"
include(tssos_dir*"tssos_model.jl")

# random example
using Random
Random.seed!(1234); # generate repeatable random numbers

d = 3; n = 5;
v1 = rand(5)
v1[4] = 0.0; v1[5] = 0.0
v2 = rand(5)
v2[1] = 0.0; v2[5] = 0.0
v3 = rand(5)
v3[2] = 0.0; v3[3] = 0.0
v4 = zeros(5); v4[1] = 3 # the last column records the order of this tensor
vmat = hcat(v1,v2,v3,v4)

A = matrix_IO.get_tms_from_factors(vmat) # the tms corresponding to the tensor

t₀ = Int(ceil( (d+1)/2 )) # the initial relaxation order
len = binomial(n+t₀,t₀)
G = rand(len,len)
F = G'*G

lvt = t₀

#----------------TSSOS model--------------------------------------
cpu_tssos = @elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F) end
mom_vals_tssos = mom_mat_tssos
MomMat = TSSOS.get_moment_matrix(mom_vals_tssos, info1)

tol = 1e-6
#sol = extract_solutions_robust(n, lvt, MomMat[1]; tol=tol)
@polyvar x[1:n]
cents_all = extract_solutions([zeros(n)'*x],x,lvt, 0, MomMat[1]; tol=tol)
# compute the weights
I₂ₜ = moments.make_mon_expo(n, d; isle=false)
mom = zeros(length(I₂ₜ), length(cents_all))
for i in 1:length(cents_all) 
    mom[:,i] = map(a->prod(cents_all[i].^a),I₂ₜ)
end
weights = mom \ [A[item] for item in moments.get_subscripts(d,n)]
A_ext_tssos_ideal = my_extract_atoms.recon_ten_id(cents_all,weights,d)
error_tssos_ideal = sum([abs(A_ext_tssos_ideal[key]-A[key]) for key in keys(A)])
println("tssos_ideal_recover_error: ", error_tssos_ideal)

#=#----------------------------ideal sparse model---------------------------------------------
cpu_sp = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"sp") end

n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(mom_mat_sp)

try
    println("----sp extract-----")
    ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-4, true)

    cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp,A,n,d)
    A_ext_sp                = extract_atoms.recon_ten_sp(cents_sp, weights_sp, d)

    error_sp = sum([abs(A_ext_sp[key]-A[key]) for key in keys(A)])
    println("sp_recover_error: ", error_sp)
catch
    println("No atom could be extracted")
end

#-- ideal dense model----------------------------------------------------------------------
cpu_id = @elapsed begin ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,F,"id") end

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