function [Q,T] = orderschur(A)
% [Q,T] = orderschur(A) computes the real Schur decomposition of 
% a matrix with real eigenvalues sorted in increasing order along
% the diagonal of T
%
% Algorithm: perform unordered Schur decomposition with Matlab's schur
% function, then order the eigenvalues with Givens rotations
% as the Algorithm 7.6.1 in [Golub, Van Loan. Matrix computations., 4nd]

[Q,T] = schur(A);
n = size(A,1);

order = 0;
while ~order % bubble sort
    order = 1;
    for k=1:n-1
        if T(k,k) > T(k+1,k+1)
            order = 0; % indicate the swap happens
            % Givens rotation
            [c,s] = givens(T(k,k+1),T(k+1,k+1)-T(k,k));
            T(k:k+1,k:n) = [c s;-s c]'*T(k:k+1,k:n);
            T(1:k+1,k:k+1) = T(1:k+1,k:k+1)*[c s;-s c];
            Q(1:n,k:k+1) = Q(1:n,k:k+1)*[c s;-s c];
        end
    end
end

    function [c,s] = givens(a,b)
    % Givens rotation for ordered Schur decomposition
    if b == 0
      c = 1; s = 0;
    else
      if abs(b) > abs(a)
        t = -a/b; s = 1/sqrt(1+t^2); c = s*t;
      else
        t = -b/a; c = 1/sqrt(1+t^2); s = c*t;
      end
    end
    end
end