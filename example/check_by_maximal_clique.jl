src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"matrix_IO.jl")
include(src_dir*"moments.jl")

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


n = parse(Int64,split(A_name,'_')[2][2:end])
d = parse(Int64,split(A_name,'_')[3][2:end])

@elapsed begin
mc = moments.get_maximal_cliques(A,n,d)
might_CP = all(a->!isempty(a),[filter(c->issubset(item,c),mc) for item in moments.get_nonzero_entries(A)])
end