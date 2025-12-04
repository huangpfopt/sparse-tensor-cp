function [rs,check] = rank_check(n,lvt,d,Mw,tol)
rs = [];
oldrank = 0;
check = false;
for j = ceil(d/2)-1:lvt
    sz = nchoosek(n+j,j);

    [~,S] = svd(Mw(1:sz,1:sz)); % [~,S] = eig(Mw(1:sz,1:sz)) at least results of eig and svd are similar when compute the rank
    S = diag(S);
    drop = find(S>=tol);
    r = length(drop);
    rs = [rs,r];
    if r == oldrank
        check = true;
        break
    else
        oldrank = r;
    end
end

end