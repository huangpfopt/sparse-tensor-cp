%% load data
data_path = "..\assets\data\literature_non\";
data_tensor_dict_noncp = data(data_path);

dmdk = data_tensor_dict_noncp.keys;

k = 2;
T_name = dmdk{k}; 
T = data_tensor_dict_noncp(dmdk{k});

%% SDP relaxation
spl = split(T_name,'_');
n = eval(spl{2}(2:end));% dimension of tensor
d = eval(spl{3}(2:end));% order of tensor
t0 = ceil((d+1)/2);
%len = nchoosek(t0+n,n);
%G = rand(len,len);
%F = G'*G;
F = readmatrix(sprintf("F\\F_%s_%d.csv",T_name(1:end-4),t0));

lvt = t0; % the relaxation order

%% gloptipoly mosek model
[output_mosek,mosek_info] = glopti_model_mosek(T,n,d,lvt,F); 
fprintf("mosek_info:\n")
disp(mosek_info)
fprintf("objective value is:%f\n", output_mosek.obj)
fprintf("sdptime is:%f\n", output_mosek.sdptime)
fprintf("cputime is:%f\n", output_mosek.cpu)

Mw_mosek = output_mosek.Mw;

%% check rank condition
tol = 1e-6;
[rs_mosek,check_mosek] = rank_check(n,lvt,d,Mw_mosek,tol);
fprintf("rs_mosek:\n")
disp(rs_mosek);
fprintf("check_mosek:\n")
disp(check_mosek);

%% extract atoms
try
    fprintf("----svd extract---\n")
    [cents_mosek, weights_mosek,error_mosek] = get_atoms(T,n,lvt,Mw_mosek,tol);
    fprintf("cents_mosek:\n")
    disp(cents_mosek);
    fprintf("weights_mosek:\n")
    disp(weights_mosek);
    fprintf("mosek recover error: %.4e\n", error_mosek);

catch
    fprintf("No atoms can be extracted\n")
end

%% gloptipoly sedumi model
[output,sedumi_info] = glopti_model(T,n,d,lvt,F); %sedumi_info.feasration < 0 indicates that infeasible (ref Ni's code)
fprintf("sedumi_info:\n")
disp(sedumi_info)
fprintf("objective value is:%f\n", output.obj)
fprintf("sdptime is:%f\n", output.sdptime)
fprintf("cputime is:%f\n", output.cpu)

Mw = output.Mw;

%% check rank condition
tol = 1e-6;
[rs,check] = rank_check(n,lvt,d,Mw,tol);
fprintf("rs:\n")
disp(rs);
fprintf("check:\n")
disp(check);

%% extract atoms
try
    fprintf("----svd extract---\n")
    [cents, weights,error] = get_atoms(T,n,lvt,Mw,tol);
    fprintf("cents:\n")
    disp(cents);
    fprintf("weights:\n")
    disp(weights);
    fprintf("recover error: %.4e\n", error);

catch
    fprintf("No atoms can be extracted\n")
end

