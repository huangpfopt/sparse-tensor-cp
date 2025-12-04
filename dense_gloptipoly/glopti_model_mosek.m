function [output,info] = glopti_model_mosek(T,n,d,lvt,F)
start = tic;
mset clear % delete all existing Gloptipoly variables from the matlab workspace
mset('yalmip',true)
mset(sdpsettings('solver','mosek'));
mpol('x',n) % define the n-dimensional variables of polynomials, alternatively, mpol x n

Amon = mmon(x,d,d); % return the vector of all monomials of degree d, e.g. x(1)^d, x(1)^(d-1)x(2),x(1)^(d-1)x(3),...,x(n)^d
conmom = [mom(Amon) == T.observes];
K = [x'*x-1==0, x'>=0];

t0 = ceil((d+1)/2);
bracx = mmon(x,0,t0); % return the vector of all monomials of degree \le k0
P = msdp(min(mom(bracx'*F*bracx)),K,conmom,lvt); %transform the linear moment optimization to a sdp relaxation problem
[F2,h,y] = myalmip(P); % return the yalmip form of P

substart = tic;
info = solvesdp(F2,h);
sdptime = toc(substart);
fprintf('\n mosek solve time for order %d is %.4f\n', lvt, sdptime);

cpu = toc(start);
ysol = double(y);
[A,b,c,S] = msedumi(P);
dvar = c-A'*ysol;
Mw = mat(dvar(S.f+1:S.f+S.s(1)^2)); % the moment metrix, mat is a command in sedumi

output.sdptime = sdptime;
output.cpu = cpu;
output.Mw = Mw;
output.obj = sum(sum(F.*Mw(1:size(F,1),1:size(F,1)) ) );
output.h = double(h);
output.ysol = ysol;
end