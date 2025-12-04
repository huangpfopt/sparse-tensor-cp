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

mat_dir = data_dir*mat_types[3]*"\\"
#A_name = "ex031_n10_d4_p0.4__.csv"
A_name = "ex030_n8_d8_p0.98__.csv" #TSSOS model will report Mosek.MosekError(1051,"")
A = matrix_IO.load_tmses(mat_dir*"\\R_"*A_name)

#----------------------------------SDP relaxation--------------------------------------
n = parse(Int64,split(A_name,"_")[2][2:end])
d = parse(Int64,split(A_name,"_")[3][2:end])
t0 = Int(ceil((d+1)/2))
len = binomial(t0+n,t0)
G = rand(len,len)
F = G'*G 


lvt = t0

#=#----------------TSSOS model--------------------------------------
println("---------TSSOS begin------------")
cpu_tssos = @elapsed begin γₜᵗˢˢᵒˢ, mom_mat_tssos, info1 = tssos_model.get_γₜᵗˢˢᵒˢ(A,n,d,lvt,F) end
mom_vals_tssos = mom_mat_tssos
MomMat = TSSOS.get_moment_matrix(mom_vals_tssos, info1)

#tol = 1e-6
#sol = extract_solutions_robust(n, lvt, MomMat[1]; tol=tol)=#

#---------------check by maximal clique
clique_cpu = @elapsed begin
mc = moments.get_maximal_cliques(A,n,d)
might_CP = all(a->!isempty(a),[filter(c->issubset(item,c),mc) for item in moments.get_nonzero_entries(A)])
end
