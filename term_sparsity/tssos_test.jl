using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using TSSOS

#----------------------- Initialization ---------------------
src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"extract_atoms.jl")
include(src_dir*"cp_model.jl")

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
subscripts = moments.get_subscripts(d, n)
Avals = [A[item] for item in subscripts]

# construt the dual problem
#=
max_{p ∈ R[x]_E }<p,y>
s.t. F-p = σ₀+∑ᵢ₌₁ⁿσᵢxᵢ+τ(x'*x-1)
=#
@polyvar x[1:n] # Define the variables of polynomials
t₀ = Int(ceil( (d+1)/2 )) # the initial relaxatoin order
len = binomial(n+t₀,t₀)
bracx = vcat([MultivariatePolynomials.monomials(x, i) for i=0:t₀]...) # return all the monomials of degree ≤ t₀
G = rand(len,len)
F = bracx'*(G'*G)*bracx # construct a sos randomly

h = [x'*x-1]
g = [x...] # the inequality constraint

lvt = t₀ # the relaxation order

cpu_tssos = @elapsed begin
model = Model(optimizer_with_attributes(Mosek.Optimizer))
@variable(model, p[1:length(A)])
nonne = F - sum(p.*MultivariatePolynomials.monomials(x, d))
# QUIET：true then no process information will be output  
# CS: correlative sparse,"MF" by default (approximately smallest chordal extension), "NC" (not performing chordal extension), false (invalidating correlative sparsity exploitation)
# TS: term sparse, "block" by default (maximal chordal extension), "signsymmetry" (sign symmetries), "MD" (approximately smallest chordal extension), false (invalidating term sparsity iterations)
# SO: sparse order
info1 = add_psatz!(model, nonne, x, g, h, lvt, QUIET=false, CS=true, cliques=[], TS="block", SO=1, Groebnerbasis=false, constrs="con1")
@objective(model, Max, sum(p.*Avals))
optimize!(model)
objv = objective_value(model)
@show objv
end

# retrieve moment matrices
# mom_vals_tssos = -dual(constraint_by_name(model, "con1")) #obtain the desired tms z
MomMat = get_moment_matrix(-dual(constraint_by_name(model, "con1")), info1)

tol = 1e-6
sol = extract_solutions_robust(n, lvt, MomMat[1]; tol=tol)

#----------------------------ideal sparse model---------------------------------------------
cpu_sp = @elapsed begin ξₜᶜᵖˢᵖ, mom_mat_sp = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,G'*G,"sp") end

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
cpu_id = @elapsed begin ξₜᶜᵖⁱᵈ, mom_mat_id = cp_model.get_ξₜᶜᵖ(A,n,d,lvt,G'*G,"id") end

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
end