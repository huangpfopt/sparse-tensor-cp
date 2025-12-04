function [rs,check] = rank_check2(n,lvt,d,Mw,tol)
rs = [];
oldrank = 0;
check = false;
for j = ceil(d/2)-1:lvt
    sz = nchoosek(n+j,j);

    [~,S] = eig(Mw(1:sz,1:sz));
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