function ia = orientation_is_adjacent(comps,c,cs)
% is c adjancent to any elements in cs
% (comps, 1, [2 3])
if comps(c) == 1
    ia = any(ismember(comps(cs),[2 8])); 
elseif comps(c) == 8
    ia = any(ismember(comps(cs),[1 7])); 
else
    ia = any(ismember(comps(cs),[comps(c)-1 comps(c)+1])); 
end

