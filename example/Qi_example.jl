src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"matrix_IO.jl")

using Test

### Qi Example2(1)
v1 = zeros(10,1)
v1[2] = 1; v1[6]=1; v1[9]=1
v2 = zeros(10,1)
v2[2] = v2[8] = v2[10] = 1
v3 = zeros(10,1)
v3[3] = v3[4] = v3[5] = 1
v4 = zeros(10,1)
v4[7] = v4[9] = v4[10] = 1
v5 = zeros(10,1)
v5[1] = v5[5] = 1
v6 = zeros(10,1)
v6[2] = v6[3] = 1
v7 = zeros(10,1)
v7[5] = v7[9] = 1
v8 = zeros(10,1)
v8[2] = 1.2599
v9 = zeros(10,1); v10 = zeros(10,1); v11 = zeros(10,1)
v12 = zeros(10,1); v13 = zeros(10,1); v14 = zeros(10,1)
v9[3]=v10[4]=v11[5]=v12[6]=v13[7]=v14[8]=1
v15 = zeros(10,1)
v15[9] = 1.2599
v16 = zeros(10,1)
v16[10] = 1.2599
v17 = zeros(10,1)
v17[1] = 3;

vmat = hcat(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17)
Qi2_1 = matrix_IO.get_tms_from_factors(vmat)

### Qi Example2(2)
v1 = zeros(10,1)
v1[1] = 1; v1[5]=1; v1[10]=1
v2 = zeros(10,1)
v2[2] = v2[3] = v2[9] = 1
v3 = zeros(10,1)
v3[2] = v3[9] = v3[10] = 1
v4 = zeros(10,1)
v4[3] = v4[4] = v4[8] = 1
v5 = zeros(10,1)
v5[3] = v5[8] = v5[9] = 1
v6 = zeros(10,1)
v6[2] = v6[8] = 1
v7 = zeros(10,1)
v7[5] = v7[8] = 1
v8 = zeros(10,1)
v8[1] = 1
v9 = zeros(10,1); v10 = zeros(10,1); v11 = zeros(10,1)
v12 = zeros(10,1); v13 = zeros(10,1); v14 = zeros(10,1)
v9[2] = 1.2599
v10[3] = 1.4422
v11[4]=v12[5]=1
v13[8]=1.2599
v14[9]=1.4422
v15 = zeros(10,1)
v15[10] = 1.2599
v17 = zeros(10,1)
v17[1] = 3;

vmat = hcat(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v17)
Qi2_2 = matrix_IO.get_tms_from_factors(vmat)

### Qi Example3(2)
v1 = zeros(10,1)
v1[1] = 1; v1[3]=1; v1[8]=1; v1[9]=1;
v2 = zeros(10,1)
v2[1] = v2[2] = v2[7] = 1
v3 = zeros(10,1)
v3[2] = v3[3] = v3[6] = 1
v4 = zeros(10,1)
v4[2] = v4[6] = v4[7] = 1
v5 = zeros(10,1)
v5[7] = v5[9] =v5[10]= 1
v6 = zeros(10,1)
v6[1] = v6[3] = 1
v7 = zeros(10,1)
v7[1] = v7[5] = 1
v8 = zeros(10,1)
v8[1] = v8[8]=1
v9 = zeros(10,1); v10 = zeros(10,1); v11 = zeros(10,1)
v12 = zeros(10,1); v13 = zeros(10,1); v14 = zeros(10,1)
v9[1]=v9[9] = 1
v10[3]=v10[8]=1
v11[3]=v11[9]=1
v12[4]=v12[9]=1
v13[8]=v13[9]=1
v14[1]=1.3161
v15 = zeros(10,1)
v15[2] = 1.3161
v16 = zeros(10,1)
v16[3] = 1.3161
v17 = zeros(10,1)
v17[6] = 1.1892;
v18 = zeros(10,1); v19=zeros(10,1);v20=zeros(10,1);
v18[7] = 1.3161; v19[8]=1.1892; v20[9]=1.3161
v21=zeros(10,1);
v21[10]=1
v22=zeros(10,1);
v22[1] = 4;

vmat = hcat(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,v21,v22)
Qi3_2 = matrix_IO.get_tms_from_factors(vmat)


### Qi Example3(3)
v1 = zeros(10,1)
v1[1] = 1; v1[5]=1; v1[6]=1; v1[8]=1;
v2 = zeros(10,1)
v2[1] = v2[5] = v2[6] = v2[9] = 1
v3 = zeros(10,1)
v3[1] = v3[5] = v3[9] = v3[10] = 1
v4 = zeros(10,1)
v4[2] = v4[5] = v4[6] =v4[9] = 1
v5 = zeros(10,1)
v5[3] = v5[4] =v5[9]= 1
v6 = zeros(10,1)
v6[3] = v6[9] = v6[10] = 1
v7 = zeros(10,1)
v7[5] = v7[8] = v7[9] = 1
v8 = zeros(10,1)
v8[1] = v8[5]=1.3161
v9 = zeros(10,1); v10 = zeros(10,1); v11 = zeros(10,1)
v12 = zeros(10,1); v13 = zeros(10,1); v14 = zeros(10,1)
v9[1]=v9[6] = 1.1892
v10[1]=v10[8]=1
v11[1]=v11[9]=1.1892
v12[1]=v12[10]=1
v13[2]=v13[5]=1
v14[2]=v14[6] = 1
v15 = zeros(10,1)
v15[2]=v15[9] = 1
v16 = zeros(10,1)
v16[5] = v16[6] = 1.3161
v17 = zeros(10,1)
v17[5] = v17[8] = 1
v18 = zeros(10,1); v19=zeros(10,1);v20=zeros(10,1);
v18[5] =v18[9] = 1.3161; v19[5]=v19[10]=1; v20[6]=v20[8]=1
v21=zeros(10,1);
v21[6]=v21[9]=1.1892
v22=zeros(10,1);
v22[9] =v22[10] =  1;
v23 = zeros(10,1); v24 = zeros(10,1); v25 = zeros(10,1); v26 = zeros(10,1);
v27 = zeros(10,1); v28 = zeros(10,1); v29 = zeros(10,1); v30 = zeros(10,1);
v23[1] = 1.5651; v24[2] = 1.1892; v25[3] = 1.1892; v26[4] =1;
v27[5] = 1.7321; v28[6] = 1.5651; v29[8] = 1.3161; v30[9] = 1.7321;
v31 = zeros(10,1); v31[10] = 1.3161
v32 = zeros(10,1); v32[1] = 4;

vmat = hcat(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31,v32)
Qi3_3 = matrix_IO.get_tms_from_factors(vmat)


assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
mat_types = ["literature", "literature_non", "random"]

mat_dir = data_dir*mat_types[1]*"\\"
data_tensor_dict_cp = matrix_IO.load_tmses(mat_dir)

dmdk = data_tensor_dict_cp.keys

@test all([data_tensor_dict_cp[dmdk[4]][item]==round(Qi2_1[item]) for item in keys(Qi2_1)])
@test all([data_tensor_dict_cp[dmdk[5]][item]==round(Qi2_2[item]) for item in keys(Qi2_2)])
@test all([data_tensor_dict_cp[dmdk[6]][item]==round(Qi3_2[item]) for item in keys(Qi3_2)])
@test all([data_tensor_dict_cp[dmdk[7]][item]==round(Qi3_3[item]) for item in keys(Qi3_3)])
