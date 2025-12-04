function subscripts = get_subscripts(d,n)
% input: n is the dimension, d is the tensor order
% output: subscripts: ordered subscripts of tensors corresponding to the order of monomials

len = nchoosek(n+d-1,d);
subscripts = zeros(len,d);

v = ones(1,d);
for row = 1:len
    subscripts(row,:) = v;
    for i = d:-1:1
        if v(i) < n
            v(i) = v(i) + 1;
            for j=i+1:d
                v(j) = v(i);
            end
            break;
        end
        % only when v(i)=n, the index at the i-1-th position will change
    end
end
end