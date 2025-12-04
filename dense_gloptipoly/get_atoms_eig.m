function [cents, weights,error] = get_atoms_eig(T,n,lvt,Mw,tol)
% Step1: Compute the Choleskey decomposition of Mw, Mw=UU^T
%[U,S] = svd(Mw); % [U,S,V] = svd(Mw), U*S*V'=Mw, the diagonal of S is ordered from largest to smallest
[U,S] = eig(Mw); % [U,S] = eig(Mw), U are right eigenvector, that is A*U=U*D, the eigenvalue in S may not be sorted yet
S = diag(S);
drop = find(S>=tol);
rankMw = length(drop);
U = U(:,drop)*diag(sqrt(S(drop)));

% Step2: reduce U to the column echelon form
[UT,basis] = rref(U',tol); % basis records the row index of the leading nonzero element in each column of the reduced column-echelon form of U
U = UT';

% Step3: extract multiplication matrix for each variable
N = cell(1,n); % N{i}=M(x_iy)

% Generate powers
vpow = [];
for k=0:lvt
    vpow = [vpow;genpow(n,k)]; % genpow(3,2): [2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2]
end
BasisPow = vpow(basis,:);%obtain alpha the power index of basis monomials

for i=1:n
    ei = zeros(1,n);
    ei(i)=1;
    row = zeros(1,rankMw);
    for j = 1:rankMw
        %find the corresponding row index in U
        row(j) = findrow(BasisPow(j,:)+ei,vpow); % findrow(sub,subs) will find the row of sub in subs
    end
    N{i} = U(row,:);
end

% random combination of multiplication matrices
coef = rand(n,1);
coef = coef/sum(coef);
M = zeros(rankMw);
for i = 1:n
    M = M+coef(i)*N{i};
end

% Step4: perform ordered schur decomposition of M
[Q,~] = orderschur(M);

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