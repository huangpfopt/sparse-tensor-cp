"""
To activate enviroment:
using Pkg
Pkg.activate(".")
Pkg.instantiate()

To fun file in vscode: 
ctrl+J >Julia
> include("src\\demo.jl"), 
it seems that mosek can use multi-threads automatically in Julia mode, TSSOS can be faster than in REPL.
"unsolved": however, from the print information, optimizer thread is 14 in both case, which is weird

To run in julia REPL:
Alt+J Alt+O
> include("src\\demo.jl")
but this will lead that TSSOS cost more time, 
and performance vary from run to run, times can range from 14s to 70s for d=4 case 
(weird, might also be influenced by threads and memory occupy?)

one can Shift+Enter to run extract atoms independently, after run demo(include extract procedure) twice
since we didn't set random seed in my_extract_atoms.jl, run multiple times will give different results

It seems that my_extrat_atoms.get_atoms_robust is the most stable 
and obtains highest recovery accuracy for ideal sparse model

To evaluete the execute time:
we use @elapsed begin xxx end
and run the code twice, since it has compile time for the first time
However, it seems that the compile time only affect the first example (of random and of literature respectively) significantly.
Alternatively, one can uncomment the warm-up code in batch-run.jl line 284, it also run each example twice
"""
#----------------------------------Initializatoin--------------------------------------
src_dir = dirname(@__FILE__)*"\\" #"...\\sp_tensor_cp\\src"
include(src_dir*"matrix_IO.jl")

#-----------------------------------Generating Data------------------------------------
include(src_dir*"generate_data.jl")


#-----------------------------------Load Data------------------------------------------
println("----load data---------")
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
mat_types = ["literature", "literature_non", "random"]

mat_dir = data_dir*mat_types[1]*"\\"
data_tensor_dict_cp = matrix_IO.load_tmses(mat_dir)

mat_dir = data_dir*mat_types[2]*"\\"
data_tensor_dict_noncp = matrix_IO.load_tmses(mat_dir)

mat_dir = data_dir*mat_types[3]*"\\"
data_tensor_dict_rand = matrix_IO.load_tmses(mat_dir)


#-----------------------------------Batch computations----------------------------------
include(src_dir*"batch_run.jl")

## for maximal clique tests
println("----batch run for maximal clique tests---------")

results_dir = assets_dir*"results\\"
!isdir(results_dir) ? mkdir(results_dir) : 0
results_random_dir = results_dir*"random\\"
!isdir(results_random_dir) ? mkdir(results_random_dir) : 0
results_subdir = results_random_dir*"clique\\"
!isdir(results_subdir) ? mkdir(results_subdir) : 0

batch_run.batch_clique_save(data_tensor_dict_rand,results_subdir)

## for comparison of ideal sparse model and tssos model
println("----batch run for comparison between ideal sparse model and tssos model-----------")

results_dir = assets_dir*"results\\"
!isdir(results_dir) ? mkdir(results_dir) : 0
results_literature_dir = results_dir*"literature\\"
!isdir(results_literature_dir) ? mkdir(results_literature_dir) : 0

t = 0 # relaxation order lvt will be t0+t
batch_run.batch_comp_sparse_tssos_save(data_tensor_dict_cp, t, results_literature_dir)

dmdk = data_tensor_dict_cp.keys
tens_dict = Dict()
for k in [1,2,4,5]
    tens_dict[dmdk[k]] = data_tensor_dict_cp[dmdk[k]]
end
batch_run.batch_comp_dense_save(tens_dict, t, results_literature_dir)

t = 0
results_literature_non_dir = results_dir*"literature_non\\"
!isdir(results_literature_non_dir) ? mkdir(results_literature_non_dir) : 0
batch_run.batch_comp_sparse_tssos_save(data_tensor_dict_noncp, t, results_literature_non_dir)

batch_run.batch_comp_dense_save(data_tensor_dict_noncp, t, results_literature_non_dir)

#-----------------------------------Post Processing---------------------------------------
batch_run.make_summary_clique(results_subdir)

batch_run.make_summary_compare(results_literature_dir*"t2\\")
batch_run.make_summary_compare(results_literature_dir*"t3\\")

batch_run.make_summary_compare(results_literature_non_dir*"t2\\")
batch_run.make_summary_compare(results_literature_non_dir*"t3\\")

#-----------------------------------Extract Atoms-----------------------------------------
moments_dir = results_literature_dir*"t2\\moments\\"
batch_run.batch_extract_atoms(data_tensor_dict_cp,moments_dir)

batch_run.batch_extract_id_atoms(data_tensor_dict_cp,moments_dir)

moments_dir = results_literature_dir*"t3\\moments\\"
batch_run.batch_extract_atoms(data_tensor_dict_cp,moments_dir) 
