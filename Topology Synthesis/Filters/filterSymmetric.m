%% filterSymmetric()

function isSym = filterSymmetric(G1,G_SET)

if G1.g == 1 %G1 is the first of its generation, no comparison can be made
    isSym = false;
    return
else
    Set_nextgen = G_SET(G1.gen,:); 
end
%Indices of these graphs
V_g = find(~cellfun(@isempty,(Set_nextgen)));

%default
isSym = false;
%check with graphs in current writing gen
for g = V_g
     G2 = Set_nextgen{g};
     tf =  checkSymmetric(G1,G2);
     if tf
         isSym = true;
         "Symmetric: ["+G1.label+"; "+G2.label+"]"
         return
     end
end