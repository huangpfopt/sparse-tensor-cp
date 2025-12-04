src = dirname(@__FILE__)*"\\"
include(src*"cp_tensors.jl")
using .cp_tensors # load a module from a package: using module, while load a module from a locally defined module, using .module

assets_dir = dirname(dirname(@__FILE__))*"\\assets\\"
!isdir(assets_dir) ? mkdir(assets_dir) : 0
data_dir = assets_dir*"\\data\\"
!isdir(data_dir) ? mkdir(data_dir) : 0

tensor_types = ["literature", "literature_non", "random"]

# generate the cp tensor in lterature
ten_dir = data_dir*tensor_types[1]*"\\"
cp_tensors.generate_lit_cp_tens(ten_dir)


# generate the noncp tensor in literature
ten_dir = data_dir*tensor_types[2]*"\\"
cp_tensors.generate_lit_non_cp_tens(ten_dir)

# generate random generated sparse symmetric tensor
ten_dir = data_dir*tensor_types[3]*"\\"
# n = 12, m=16 the generate time is apprarently longer; n=12，m=26，we terminate it since the runtime is too long to stop
# generate_random_cp_tens((n,m),p)
cp_tensors.generate_random_cp_tens(ten_dir) #dimension of tensors n = [6,8,10,12,14] ; order of tensors m = [4,6,8] ; the ratio of nonzeros p = [0.4,0.6,0.8,0.9,0.98]
