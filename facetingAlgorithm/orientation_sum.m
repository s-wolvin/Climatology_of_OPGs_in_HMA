function fs = orientation_sum(comps,fracs,c)

if comps(c) == 1
    fs = max(fracs(comps==1)+fracs(comps==8),fracs(comps==1)+fracs(comps==2)); 

elseif comps(c) == 2
    fs = max(fracs(comps==2)+fracs(comps==1),fracs(comps==2)+fracs(comps==3));

elseif comps(c) == 3
    fs = max(fracs(comps==3)+fracs(comps==2),fracs(comps==3)+fracs(comps==4));

elseif comps(c) == 4
    fs = max(fracs(comps==4)+fracs(comps==3),fracs(comps==4)+fracs(comps==5));

elseif comps(c) == 5
    fs = max(fracs(comps==5)+fracs(comps==4),fracs(comps==5)+fracs(comps==6));

elseif comps(c) == 6
    fs = max(fracs(comps==6)+fracs(comps==5),fracs(comps==6)+fracs(comps==7));

elseif comps(c) == 7
    fs = max(fracs(comps==7)+fracs(comps==6),fracs(comps==7)+fracs(comps==8));

elseif comps(c) == 8 
    fs = max(fracs(comps==8)+fracs(comps==7),fracs(comps==8)+fracs(comps==1)); 
end
