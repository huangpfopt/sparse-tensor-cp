%% load data
data_path = "..\assets\data\literature\";
data_tensor_dict_cp = data(data_path);

%data_path = "..\assets\data\literature_non\";
%data_tensor_dict_noncp = data(data_path);


dmdk = data_tensor_dict_cp.keys;
%dmdk = data_tensor_dict_noncp.keys;

for k=1:length(dmdk)
%k = 1;
T_name = dmdk{k}; %'ex01_n11_d3_17Fan-Exa4.1_.csv'
T = data_tensor_dict_cp(dmdk{k});
%T = data_tensor_dict_noncp(dmdk{k}); % T.subscripts(including all subscripts by order), T.observes


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
save_path = sprintf(".\\results\\mosek_glopti_%s_%d.txt",T_name(1:end-4),lvt);
%save_path = sprintf(".\\results\\mosek_glopti_noncp%s_%d.txt",T_name(1:end-4),lvt);
diary(save_path);
diary on;

%% gloptipoly model
[output,mosek_info] = glopti_model_mosek(T,n,d,lvt,F); 
fprintf("mosek_info:\n")
disp(mosek_info)
fprintf("objective value is:%f\n", output.h)
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
    [cents, weights,error] = get_atoms(T,n,lvt,Mw,tol);
    fprintf("cents:\n")
    disp(cents);
    fprintf("weights:\n")
    disp(weights);
    fprintf("recover error: %.4e\n", error);

catch
    fprintf("No atoms can be extracted\n")
end
diary off;
end