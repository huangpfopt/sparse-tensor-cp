function index = findrow(sub,subs)
for i = 1:length(subs)
    if sub == subs(i,:)
        index = i;
        break
    end
end
end