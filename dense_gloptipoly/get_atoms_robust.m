function [cents, weights,error] = get_atoms_robust(T,n,lvt,Mw,tol)
% Generate powers
vpow = [];
for k=0:lvt
    vpow = [vpow;genpow(n,k)]; % genpow(3,2): [2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2]
end
ls = nchoosek(n+lvt-1,n);

for i=1:n
    ei = zeros(1,n);
    ei(i)=1;
    N{i} = zeros(ls);
    for j = 1:ls
        for k = 1:ls
            loc = findrow(vpow(k,:)+ei,vpow);
            N{i}(j,k) = Mw(j,loc);
        end
    end
    N{i} = triu(N{i})+triu(N{i},1)'; % triu(N{i},k), k=0包括主对角线的上三角矩阵，k=1,主对角线取0的上三角部分矩阵
end
[U,S,V] = svd(Mw(1:ls,1:ls));
S = diag(S);
drop = find(S>=tol);
rankMw = length(drop);
S = sqrt(S(drop)).^(-1);
for i = 1:n
    N{i} = diag(S)*V(:,drop)'*N{i}*U(:,drop)*diag(S);
end
coef = rand(n,1);
coef = coef/sum(coef);
M = zeros(rankMw);
for i = 1:n
    M = M+coef(i)*N{i};
end

[Q,~] = schur(M); % QTQ'=M

% Step5: rectrieve atoms
cents = zeros(n,rankMw);
for j=1:rankMw
    for i=1:n
        cents(i,j) = Q(:,j)'*N{i}*Q(:,j);
    end
end

% compute the weights
T_ext = zeros(length(T.observes),rankMw);
subscripts = T.subscripts;
for i = 1:length(T.observes)
    for j = 1:rankMw
        x = cents(:,j);
        T_ext(i,j) = prod(x(subscripts(i,:)));
    end
end
weights = T_ext\T.observes;
error = sum(abs(T.observes - T_ext*weights));

end