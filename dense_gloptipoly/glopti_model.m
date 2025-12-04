function [output,sedumi_info] = glopti_model(T,n,d,lvt,F)
start = tic;
mset clear % delete all existing Gloptipoly variables from the matlab workspace
mpol('x',n) % define the n-dimensional variables of polynomials, alternatively, mpol x n

Amon = mmon(x,d,d); % return the vector of all monomials of degree d, e.g. x(1)^d, x(1)^(d-1)x(2),x(1)^(d-1)x(3),...,x(n)^d
conmom = [mom(Amon) == T.observes];
K = [x'*x-1==0, x'>=0];

t0 = ceil((d+1)/2);
bracx = mmon(x,0,t0); % return the vector of all monomials of degree \le k0
P = msdp(min(mom(bracx'*F*bracx)),K,conmom,lvt); %transform the linear moment optimization to a sdp relaxation problem
[A,b,c,S] = msedumi(P); % return the sedumi form of P, it satisfies the dual problem form in sedumi

substart = tic;
[~,ysol,sedumi_info] = sedumi(A,b,c,S); %solve the sedumi form problem, ysol is the tms?
sdptime = toc(substart);
fprintf('\n sedumi solve time for order %d is %.4f\n', lvt, sdptime);
if sedumi_info.feasratio<0
    fprintf('\n infeasible for relaxation order %d\n', lvt);
end
cpu = toc(start);

dvar = c-A'*ysol;
Mw = mat(dvar(S.f+1:S.f+S.s(1)^2)); % the moment metrix, mat is a command in sedumi

output.sdptime = sdptime;
output.cpu = cpu;
output.Mw = Mw;
%output.obj = sum(sum(F.*Mw)); 
output.obj = sum(sum(F.*Mw(1:size(F,1),1:size(F,1)) ) );
output.ysol = ysol;
end